%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Redbird - Mesh-mode adjoint Jacobian via mmclab
%
% Demonstrates the 6th return value of rbrunforward:
%
%   [detphi, phi, ~, ~, ~, Jext] = rbrunforward(cfg);
%
% Set when cfg.nphoton is present and nargout > 5. rbrunforward then
%   - auto-fills cfg.detdir via rbgetdetdir when absent,
%   - sets cfg.outputtype='adjoint_mua_d', cfg.basisorder=1, cfg.srcid=-1,
%   - launches one batched mmclab call per wavelength packing all Ns
%     forward + Nd detector-adjoint slots into one GPU run,
%   - returns Jext = struct('mua', J_mua, 'dcoeff', J_D) in the
%     [Nsd x Nn] convention used by rbjac.
%
% Jext.mua/dcoeff are real for CW and complex for RF (cfg.omega > 0,
% supported since mmc's BLBadouel mesh kernel gained complex-weight
% tracking).
%
% This demo plots the J_mua banana profile for each source-detector pair
% and cross-checks it against the FEM rb_femjacobian formula
% (rbjac) on the same mesh.
%
% Companion demos:
%   demo_redbird_forward_mc.m - forward simulation via mmclab
%   demo_redbird_recon_mc.m   - CW DOT reconstruction via MC
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
%%   homogeneous 60 x 60 x 30 mm slab
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear cfg;

[cfg.node, cfg.face, cfg.elem] = meshabox([0 0 0], [60 60 30], 4);
cfg.seg = ones(size(cfg.elem, 1), 1);

cfg.srcpos = [30 30 0];
cfg.srcdir = [0  0  1];

% three detectors at increasing source-detector separation along +x
% (rbrunforward adds the disk-source radius cfg.detpos(:,4) = avgsize
%  automatically when routing to mmclab)
cfg.detpos = [20 30 0
              40 30 0
              45 30 0];
cfg.detdir = [0 0 1; 0 0 1; 0 0 1];

cfg.prop  = [0 0 1 1; 0.005 1 0 1.37];
cfg.omega = 0;

cfg = rbmeshprep(cfg);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   MC adjoint Jacobian via rbrunforward (6th output)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cfg.nphoton = 3e7;

fprintf('mmclab mesh adjoint (%.0e photons) ...\n', cfg.nphoton);
tic;
[~, phi_mc, ~, ~, ~, Jext] = rbrunforward(cfg);
fprintf('  done in %.2f s   Jext.mua shape: %s   Jext.dcoeff shape: %s\n', ...
        toc, mat2str(size(Jext.mua)), mat2str(size(Jext.dcoeff)));

% Jext.mua / Jext.dcoeff are [Nsd x Nn] (Nsd = Nsrc * Ndet).
% Each row is one source-detector pair's nodal Jacobian.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   FEM reference: rb_femjacobian on the same mesh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cfg_fem = rmfield(cfg, 'nphoton');
[~, phi_fem] = rbrun(cfg_fem);

sd = rbsdmap(cfg_fem);
Jmua_node_fem = rbjac(sd, phi_fem, cfg_fem.deldotdel, cfg_fem.elem, cfg_fem.evol);

fprintf('FEM Jmua_node shape: %s\n', mat2str(size(Jmua_node_fem)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   plot banana profiles, MC vs FEM, on the y=30 mid-plane
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nsd = size(Jext.mua, 1);

figure('Name', 'Mesh-mode adjoint Jacobian: MC vs FEM');
for k = 1:Nsd
    Jmc  = abs(Jext.mua(k, :))';
    Jfem = abs(Jmua_node_fem(k, :))';

    subplot(2, Nsd, k);
    plotmesh([cfg.node, log10(Jmc + 1e-14)], cfg.elem, 'y=30', ...
             'facecolor', 'interp', 'linestyle', 'none');
    view([0 1 0]);
    colorbar;
    title(sprintf('MC  |J_{\\mu_a}|  pair %d', k));

    subplot(2, Nsd, Nsd + k);
    plotmesh([cfg.node, log10(Jfem + 1e-14)], cfg.elem, 'y=30', ...
             'facecolor', 'interp', 'linestyle', 'none');
    view([0 1 0]);
    colorbar;
    title(sprintf('FEM |J_{\\mu_a}|  pair %d', k));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   quantitative agreement: log10 correlation per pair
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\nMC vs FEM J_mua agreement (log10-log10 corr coef):\n');
for k = 1:Nsd
    mask = (abs(Jext.mua(k, :)) > 1e-12) & (abs(Jmua_node_fem(k, :)) > 1e-12);
    if any(mask)
        cc = corrcoef(log10(abs(Jext.mua(k, mask))), ...
                      log10(abs(Jmua_node_fem(k, mask))));
        fprintf('  pair %d:  %.4f  (expected close to 1.0)\n', k, cc(1, 2));
    end
end
