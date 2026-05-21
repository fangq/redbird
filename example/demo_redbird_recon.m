%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Redbird - A Diffusion Solver for Diffuse Optical Tomography,
%      Copyright Qianqina Fang, 2018
%
% Continuous-Wave (CW) reconstruction of absorption (mua) target
% (streamlined version by calling rbrunrecon)
%
% This file is part of Redbird URL:http://mcx.sf.net/mmc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (~exist('rbrun', 'file'))
    addpath(fullfile(pwd, '../matlab'));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   prepare simulation input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear cfg cfg0 recon;

s0 = [70, 50, 20];

[nobbx, fcbbx] = meshabox([40 0 0], [160, 120, 60], 10);
[nosp, fcsp] = meshasphere(s0, 5, 1);
[no, fc] = mergemesh(nobbx, fcbbx, nosp, fcsp);

[cfg0.node, cfg0.elem] = s2m(no, fc(:, 1:3), 1, 40, 'tetgen', [41 1 1; s0]);
nn = size(cfg0.node, 1);
cfg0.seg = cfg0.elem(:, 5);
cfg0.srcdir = [0 0 1];

[xi, yi] = meshgrid(60:20:140, 20:20:100);
cfg0.srcpos = [xi(:), yi(:), zeros(numel(yi), 1)];
cfg0.detpos = [xi(:), yi(:), 60 * ones(numel(yi), 1)];
cfg0.detdir = [0 0 -1];

cfg0.prop = [
             0 0 1 1
             0.008 1 0 1.37
             0.016 1 0 1.37
            ];

cfg0.omega = 0;

cfg = cfg0;

cfg0 = rbmeshprep(cfg0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Run forward for the heterogeneous domain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

detphi0 = rbrun(cfg0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Reset the domain to a homogeneous medium for recon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create forward mesh for reconstruction - if it has a different density
% compared to the forward mesh used for generating data, it may introduce
% numerical error - so, set density to 10 gives the best result, need
% to debug this further

[node, face, elem] = meshabox([40 0 0], [160, 120, 60], 10);
cfg = rbsetmesh(cfg, node, elem, cfg.prop, ones(size(node, 1), 1));

sd = rbsdmap(cfg);

% create coarse reconstruction mesh
[recon.node, face, recon.elem] = meshabox([40 0 0], [160, 120, 60], 20);
[recon.mapid, recon.mapweight] = tsearchn(recon.node, recon.elem, cfg.node);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Stage 1: bulk mua/musp fit (single segment)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Start from a deliberately poor initial guess. rbrun mode='bulk' sets
% recon.seg = ones(size(recon.node)), collapsing the Jacobian columns into
% a single per-segment unknown via rbmasksum. Each Gauss-Newton step then
% adjusts one global mua (and musp), giving a fast sanity check on the
% FEM forward + adjoint pipeline before per-node imaging starts.

recon.bulk = struct('mua', 0.003, 'musp', 0.6);

fprintf('=== Stage 1: bulk mua/musp fit ===\n');
newrecon = rbrun(cfg, recon, detphi0, sd, 'mode', 'bulk', 'lambda', 1e-3, 'maxiter', 5);
fprintf('Bulk fit result: mua = %g  musp = %g\n', newrecon.prop(2, 1), newrecon.prop(2, 2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Stage 2: per-node image recon seeded with stage-1 bulk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Re-seed recon.bulk from the stage-1 fit and let rbrun mode='image' rebuild
% per-node recon.prop (size(recon.node) x 4) and cfg.prop (size(cfg.node) x 4)
% from those bulk values.

recon.bulk = struct('mua', newrecon.prop(2, 1), 'musp', newrecon.prop(2, 2));
if (isfield(recon, 'prop'))
    recon = rmfield(recon, 'prop');   % let rbrun rebuild per-node from new bulk
end

fprintf('=== Stage 2: per-node image reconstruction ===\n');
[newrecon, resid, newcfg] = rbrun(cfg, recon, detphi0, sd, 'mode', 'image', 'mex', 0, 'lambda', 1e-4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Plotting results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plotmesh([newcfg.node, newcfg.prop(:, 1)], newcfg.elem, 'z=20', 'facecolor', 'interp', 'linestyle', 'none');
hold on;
plotmesh([newcfg.node, newcfg.prop(:, 1)], newcfg.elem, 'x=70', 'facecolor', 'interp', 'linestyle', 'none');
view(3);
