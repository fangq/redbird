%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Redbird - Voxel-grid Monte Carlo (mcxlab) adjoint + matrix-free LSQR
%
% Demonstrates the mcxlab-grid Jacobian path added to rbrunforward:
%
%   [detphi, phi, ~, ~, ~, Jext] = rbrunforward(cfg);
%
% engaged when cfg.nphoton AND cfg.vol are set together.  rbrunforward
% auto-fills cfg.detdir via rbgetdetdir_vol, sets cfg.srcid = -1 and
% cfg.outputtype = 'adjoint_mua_d', and returns Jext = struct('mua', ...
% 'dcoeff', ...) with .mua of shape (Nx, Ny, Nz, Ns*Nd).
%
% That 4D Jacobian is too large for the normal-equation form
% (J' * J would be Nv x Nv with Nv ~ 1e6).  rbreglsqr wraps it via
% rbjacop into a matvec/rmatvec operator and solves
%     min || J * delta_mu - (y_meas - y_model) ||_2
% via MATLAB's built-in LSQR.  Early stopping acts as regularization;
% no lambda to tune.
%
% This demo performs ONE Gauss-Newton step (linearization around the
% background mua) to keep the runtime small.  For a full iterative
% reconstruction loop, wire cfg.muavol back into the next mcxlab call's
% cfg.vol (per-voxel property encoding) and re-run rbrunforward.
%
% Requires: mcxlab on the MATLAB path; an mcx tree that supports
% outputtype='adjoint_mua_d' and srcid=-1 / -2.
%
% Companion demos:
%   demo_redbird_forward_mc.m / demo_redbird_jacobian_mc.m for the
%   mmclab (tet mesh) analogue.
%
% This file is part of Redbird URL:http://mcx.space/#redbird
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (~exist('rbrun', 'file'))
    addpath(fullfile(pwd, '../matlab'));
end
if (~exist('mcxlab', 'file'))
    error('mcxlab is not on the MATLAB path; this demo requires MCX');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   60 x 60 x 30 voxel slab (1 mm voxels)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear cfg recon;

Nx = 60;
Ny = 60;
Nz = 30;
cfg.vol  = uint8(ones(Nx, Ny, Nz));      % single tissue label
cfg.unitinmm = 1.0;

% three sources along the +x face of the slab (z = 0)
cfg.srcpos = [20 30 0; 30 30 0; 40 30 0];
cfg.srcdir = [0  0  1];
cfg.srctype = 'pencil';

% four detectors at increasing s-d separation, same surface
cfg.detpos = [15 30 0 1.5
              25 30 0 1.5
              35 30 0 1.5
              45 30 0 1.5];

cfg.prop = [0 0 1 1
            0.005 1 0 1.37];

cfg.tstart  = 0;
cfg.tend    = 5e-9;
cfg.tstep   = 5e-9;     % single-gate CW
cfg.gpuid   = 1;
cfg.autopilot = 1;
cfg.maxdetphoton = 1e6;
cfg.isreflect = 1;

Nsrc = size(cfg.srcpos, 1);
Ndet = size(cfg.detpos, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Baseline forward + adjoint Jacobian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cfg.nphoton = 5e7;

fprintf('mcxlab adjoint (%.0e photons) ...\n', cfg.nphoton);
tic;
[detphi, phi, ~, ~, ~, Jext] = rbrunforward(cfg);
fprintf('  done in %.2f s\n', toc);

fprintf('  detphi shape : %s   (Ndet=%d x Nsrc=%d)\n', mat2str(size(detphi)), Ndet, Nsrc);
fprintf('  phi shape    : %s   (Nx,Ny,Nz,Ns+Nd)\n',     mat2str(size(phi)));
fprintf('  Jext.mua shp : %s   (Nx,Ny,Nz,Ns*Nd)\n',     mat2str(size(Jext.mua)));
if (isfield(Jext, 'dcoeff'))
    fprintf('  Jext.dcoeff  : %s\n',                        mat2str(size(Jext.dcoeff)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Synthesize a "measurement" by perturbing one voxel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Add a small mua bump at voxel (30, 30, 15) of magnitude 0.005, then
% pretend the linearized model y_meas = y_model + J * delta_mu_true is
% the truth.  This is the cleanest possible test: the LSQR solve
% should recover delta_mu_true exactly (modulo MC noise + early-stop
% truncation).

delta_mu_true = zeros(Nx, Ny, Nz);
delta_mu_true(28:32, 28:32, 14:16) = 0.005;

% J has shape (Nx, Ny, Nz, Ns*Nd); reshape to 2D for one matvec
J2 = reshape(double(Jext.mua), Nx * Ny * Nz, Nsrc * Ndet);
y_meas = detphi(:) + J2.' * delta_mu_true(:);   % linearized "measurement"

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Inverse step via rbreglsqr (matrix-free LSQR)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r = y_meas - detphi(:);                  % residual (Ns*Nd vector)

fprintf('\nrbreglsqr ...\n');
tic;
[delta_mu_rec, info] = rbreglsqr(Jext.mua, r, 'maxit', 100, 'tol', 1e-8);
fprintf('  done in %.2f s   LSQR iters: %d   relres: %.3e   adjoint err: %.2e\n', ...
        toc, info.iter, info.relres, info.adjoint_err);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Compare reconstruction to ground truth
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Locate ground-truth and recovered bumps on the z = 15 mid-slice
true_slice = delta_mu_true(:, :, 15);
rec_slice  = delta_mu_rec(:, :, 15);

% energy preservation (sum of recon should be close to sum of truth)
sum_true = sum(delta_mu_true(:));
sum_rec  = sum(delta_mu_rec(:));
peak_true = max(true_slice(:));
peak_rec  = max(rec_slice(:));

fprintf('\nGround-truth sum: %.3e   Reconstructed sum: %.3e   ratio: %.3f\n', ...
        sum_true, sum_rec, sum_rec / sum_true);
fprintf('Ground-truth peak: %.3e   Recon peak (z=15 slice): %.3e\n', peak_true, peak_rec);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Plot truth vs reconstruction on the z = 15 slice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Name', 'mcxlab voxel-grid mua reconstruction (1 Gauss-Newton step)');
subplot(1, 2, 1);
imagesc(true_slice');
axis image;
colorbar;
set(gca, 'YDir', 'normal');
title('ground-truth delta mua  (z=15)');
xlabel('x (voxel)');
ylabel('y (voxel)');

subplot(1, 2, 2);
imagesc(rec_slice');
axis image;
colorbar;
set(gca, 'YDir', 'normal');
title('recovered delta mua via LSQR  (z=15)');
xlabel('x (voxel)');
ylabel('y (voxel)');
