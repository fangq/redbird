%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Redbird - Monte Carlo (mmclab) based CW DOT reconstruction
%
% Continuous-Wave (CW) reconstruction of an absorption (mua) inclusion
% using mmc as the forward-and-adjoint engine. Mirrors the streamlined
% FEM example demo_redbird_recon.m but adds cfg.nphoton to route through
% mmclab. The mesh-mode adjoint Jacobian returned by mmclab is consumed
% by rbrunrecon directly, bypassing rbjac.
%
% Requires: mmclab on the MATLAB path, built against an mmc tree that
% includes PR#1+PR#2+PR#3 (mesh adjoint + per-node prop). Photon counts
% of order 1e7 keep the per-iteration Monte Carlo noise small enough that
% Gauss-Newton converges in a handful of steps.
%
% This file is part of Redbird URL:http://mcx.space/#redbird
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (~exist('rbrun', 'file'))
    addpath(fullfile(pwd, '../matlab'));
end
if (~exist('mmclab', 'file'))
    error('mmclab is not on the MATLAB path; this demo requires MMC (mmc/mmclab)');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   ground-truth heterogeneous domain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear cfg cfg0 recon;

s0 = [70, 50, 20];

[nobbx, fcbbx] = meshabox([40 0 0], [160, 120, 60], 10);
[nosp, fcsp]   = meshasphere(s0, 5, 1);
[no, fc]       = mergemesh(nobbx, fcbbx, nosp, fcsp);

[cfg0.node, cfg0.elem] = s2m(no, fc(:, 1:3), 1, 40, 'tetgen', [41 1 1; s0]);
cfg0.seg    = cfg0.elem(:, 5);
cfg0.srcdir = [0 0 1];

[xi, yi] = meshgrid(60:20:140, 20:20:100);
cfg0.srcpos = [xi(:), yi(:), zeros(numel(yi), 1)];
cfg0.detpos = [xi(:), yi(:), 60 * ones(numel(yi), 1)];
cfg0.detdir = [0 0 -1];

cfg0.prop = [
             0     0   1   1
             0.008 1   0   1.37   % background
             0.016 1   0   1.37   % absorbing inclusion (2x background mua)
            ];

cfg0.omega = 0;

cfg = cfg0;
cfg0 = rbmeshprep(cfg0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Synthesize measurement data with the FEM forward
%%   (a "perfect" measurement; in practice this would be
%%   real boundary-detector readings or mmclab forward)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

detphi0 = rbrun(cfg0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Forward mesh + coarse recon mesh for the inversion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[node, face, elem] = meshabox([40 0 0], [160, 120, 60], 10);
cfg = rbsetmesh(cfg, node, elem, cfg.prop, ones(size(node, 1), 1));

% MC path: route cfg.nphoton through to mmclab. detdir is auto-filled
% from the surface mesh by rbrunforward via rbgetdetdir if absent.
cfg.nphoton = 1e7;
cfg = rmfield(cfg, 'detdir');   % let rbgetdetdir compute inward normals

sd = rbsdmap(cfg);

% coarse reconstruction mesh
[recon.node, face, recon.elem] = meshabox([40 0 0], [160, 120, 60], 20);
[recon.mapid, recon.mapweight]  = tsearchn(recon.node, recon.elem, cfg.node);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Streamlined MC-based reconstruction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize to homogeneous background
recon.prop = cfg.prop(ones(size(recon.node, 1), 1) + 1, :);
cfg.prop   = cfg.prop(ones(size(cfg.node, 1), 1) + 1, :);   % per-node prop
cfg        = rmfield(cfg, 'seg');

% Per-node cfg.prop is auto-detected by rbrunforward and split into
% cfg.nodemua (and cfg.nodemusp for RF) for the mmc kernel's per-node
% global-memory path. Reconstruction loop is unchanged from FEM mode -
% rbrunrecon picks up the Jacobian from rbrunforward's 6th output (Jext)
% and skips the rbjac call.

[newrecon, resid, newcfg] = rbrun(cfg, recon, detphi0, sd, 'mex', 0, 'lambda', 1e-4, 'maxiter', 5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Plot reconstructed mua slices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
plotmesh([newcfg.node, newcfg.prop(:, 1)], newcfg.elem, 'z=20', 'facecolor', 'interp', 'linestyle', 'none');
hold on;
plotmesh([newcfg.node, newcfg.prop(:, 1)], newcfg.elem, 'x=70', 'facecolor', 'interp', 'linestyle', 'none');
view(3);
title('MC-based mua reconstruction (5 Gauss-Newton iterations, 1e7 photons/iter)');
colorbar;
