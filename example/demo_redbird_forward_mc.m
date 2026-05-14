%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Redbird - Forward simulation via mmclab (Monte Carlo)
%
% rbrunforward routes to mmclab automatically when cfg.nphoton is set
% (mmclab branch); without cfg.nphoton it uses the diffusion FEM
% solver. This demo runs both on the same tet mesh and compares the
% boundary fluence at a set of detectors.
%
% The MC branch passes cfg.srcpos as an Nsrc-by-3 matrix in a single
% batched mmclab call. flux.data is reshaped back into the redbird
% [Nnode x (Nsrc+Ndet)] convention (forward slots followed by
% detector-as-adjoint slots when cfg.detdir is set).
%
% Companion demos:
%   demo_redbird_jacobian_mc.m  - mesh-mode adjoint Jacobian
%   demo_redbird_recon_mc.m     - CW DOT reconstruction via MC
%
% Requires: mmclab on the MATLAB path.
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
%%   shared geometry: 60 x 60 x 30 mm box
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear cfg cfg_mc cfg_fem;

[cfg.node, cfg.face, cfg.elem] = meshabox([0 0 0], [60 60 30], 4);
cfg.seg    = ones(size(cfg.elem, 1), 1);

% three source positions spread along the +x face of the slab
cfg.srcpos = [20 30 0; 30 30 0; 40 30 0];
cfg.srcdir = [0  0  1];                  % broadcast to all three

% detectors at increasing s-d separation on the SAME (z=0) surface
cfg.detpos = [15 30 0; 25 30 0; 35 30 0; 45 30 0];
cfg.detdir = [0 0 1; 0 0 1; 0 0 1; 0 0 1];     % inward normals (for MC adjoint)

cfg.prop  = [0 0 1 1; 0.005 1 0 1.37];
cfg.omega = 0;                            % CW

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   FEM forward (no cfg.nphoton)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cfg_fem = rbmeshprep(cfg);

fprintf('FEM forward ...\n');
tic;
[detphi_fem, phi_fem] = rbrun(cfg_fem);
fprintf('  done in %.2f s   detphi shape: %s\n', toc, mat2str(size(detphi_fem)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   MC forward via mmclab (cfg.nphoton triggers the MC branch)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cfg_mc          = cfg;
cfg_mc.nphoton  = 1e7;
% time grid: default is single-gate CW; rbrunforward fills tstart/tend/tstep

fprintf('mmclab forward (1e7 photons) ...\n');
tic;
[detphi_mc, phi_mc] = rbrunforward(cfg_mc);
fprintf('  done in %.2f s   detphi shape: %s   phi shape: %s\n', ...
        toc, mat2str(size(detphi_mc)), mat2str(size(phi_mc)));

% phi_mc is [Nnode x (Nsrc+Ndet)] -- forward fluences in columns 1..Nsrc,
% detector-as-adjoint fluences in columns Nsrc+1..end (cfg.detdir present).
Nsrc = size(cfg.srcpos, 1);
Ndet = size(cfg.detpos, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   detector-fluence agreement
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Reshape both to [Ndet x Nsrc] for direct comparison.
det_fem = reshape(detphi_fem, Ndet, Nsrc);
det_mc  = detphi_mc;        % already [Ndet x Nsrc] from rbrunforward MC branch

fprintf('\nDetector fluence agreement (MC vs FEM, log10):\n');
for s = 1:Nsrc
    for d = 1:Ndet
        if (det_fem(d, s) > 0 && det_mc(d, s) > 0)
            ratio = det_mc(d, s) / det_fem(d, s);
            fprintf('  S%d-D%d  FEM %.3e   MC %.3e   ratio %.3f\n', ...
                    s, d, det_fem(d, s), det_mc(d, s), ratio);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   per-source fluence cross-section on y=30 mid-plane
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Name', 'rbrunforward: FEM vs MC, per-source fluence');
for s = 1:Nsrc
    subplot(2, Nsrc, s);
    plotmesh([cfg.node, log10(abs(full(phi_fem(:, s))) + 1e-12)], cfg.elem, ...
             'y=30', 'facecolor', 'interp', 'linestyle', 'none');
    view([0 1 0]);
    colorbar;
    title(sprintf('FEM phi  S%d at (%g, %g, %g)', s, cfg.srcpos(s, :)));

    subplot(2, Nsrc, Nsrc + s);
    plotmesh([cfg.node, log10(abs(phi_mc(:, s)) + 1e-12)], cfg.elem, ...
             'y=30', 'facecolor', 'interp', 'linestyle', 'none');
    view([0 1 0]);
    colorbar;
    title(sprintf('MC phi   S%d', s));
end
