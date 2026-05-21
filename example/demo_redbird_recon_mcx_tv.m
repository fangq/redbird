%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Redbird - TV-regularized iterative voxel-grid mua recon via mcxlab
%
% Same outer Gauss-Newton structure as demo_redbird_recon_mcx_iter.m,
% but the inner J*delta_mu = r least-squares step is solved with
% rbregtv (TV-regularized FISTA) instead of rbreglsqr (minimum-L2-norm
% LSQR).
%
% Motivation:  with a sparse DOT optode array (Nsd = 12 here, Nv = 1.08e5
% voxels) the system is grossly underdetermined.  LSQR's min-L2-norm
% solution lives in the 12-dim row span of J and is forced to a single
% smooth blob even when the truth has multiple compact inclusions.  TV
% regularization adds a piecewise-constant prior that pulls the recon
% toward sparse-gradient (i.e., blocky) solutions, which is a much better
% match to the underlying anatomy.
%
% Tuning note:  alpha is the only knob.  Start at ~1e-4 .. 1e-3 relative
% to ||J^T r||_inf; raise alpha if the recon is still a single smooth
% blob, lower it if the recon stays at zero or develops staircase
% artifacts inside individual inclusions.  An L-curve sweep is the
% rigorous way -- see Hansen & O'Leary 1993.
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
%%   geometry: 60 x 60 x 30 voxel slab
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear cfg cfg_truth;

Nx = 60;
Ny = 60;
Nz = 30;

mua_bg = 0.005;
musp_bg = 1.0;
n_bg = 1.37;

% ground-truth perturbation: two inclusions (same as iter demo)
delta_mu_true = zeros(Nx, Ny, Nz);
delta_mu_true(10:15, 28:32, 14:16) = 0.005;
delta_mu_true(28:32, 38:42, 14:16) = 0.008;

% sources / detectors: 2D grid on z=0 covering both inclusions
%   - the original 1-D y=30 line array could only "see" one of the two
%     truth inclusions: the y=38..42 blob lay outside the banana-shaped
%     sensitivity region of any y=30 source-detector pair, so its
%     Jacobian columns were ~0 and no regularizer (LSQR or TV) could
%     recover it.
%   - the widened 4 x 4 source grid + 5 x 5 detector grid (Nsd = 400)
%     covers x in [5, 55], y in [5, 55] so both inclusions, including
%     truth 1 near the left edge (x = 10..15), sit well inside the
%     banana coverage of multiple src-det pairs.
[Sx, Sy] = ndgrid([10 25 40 55], [10 25 40 55]);
cfg_template.srcpos = [Sx(:), Sy(:), zeros(numel(Sx), 1)];
cfg_template.srcdir = [0  0  1];
cfg_template.srctype = 'pencil';
[Dx, Dy] = ndgrid([5 17 30 42 55], [5 17 30 42 55]);
cfg_template.detpos = [Dx(:), Dy(:), zeros(numel(Dx), 1), 1.5 * ones(numel(Dx), 1)];
cfg_template.prop = [0 0 1 1
                     mua_bg musp_bg 0 n_bg];
cfg_template.unitinmm = 1.0;
cfg_template.tstart = 0;
cfg_template.tend = 5e-9;
cfg_template.tstep = 5e-9;
cfg_template.gpuid = 1;
cfg_template.autopilot = 1;
cfg_template.maxdetphoton = 1e6;
cfg_template.isreflect = 1;

Nsrc = size(cfg_template.srcpos, 1);
Ndet = size(cfg_template.detpos, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Synthesize "measurements" via the full nonlinear forward
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mua_truth = single(mua_bg + delta_mu_true);
cfg_truth = cfg_template;
cfg_truth.vol = reshape(mua_truth, 1, Nx, Ny, Nz);
cfg_truth.nphoton = 5e7;

fprintf('Synthesizing target measurements (mcxlab, true mua) ...\n');
tic;
[y_meas, ~] = rbrunforward(cfg_truth);
fprintf('  done in %.2f s   y_meas shape: %s\n', toc, mat2str(size(y_meas)));

y_meas = y_meas(:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Gauss-Newton loop with TV-regularized inner solve
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mua_est = single(mua_bg * ones(Nx, Ny, Nz));
cfg = cfg_template;
cfg.nphoton = 5e7;

maxiter = 4;
resid_history = zeros(maxiter, 1);

% TV regularization weight; tune on the first iter's ||J^T r||_inf scale
% (positivity constraint takes over the spurious-negative role that TV
% used to play, so we can afford a much smaller alpha_rel here)
alpha_rel = 1e-10;

% best-iterate tracking: GN without a line search can diverge once the
% residual approaches the MC noise floor (a too-aggressive delta_mu in
% the linearized solve pushes mua into a regime where the linearization
% breaks down).  We keep mua_est at the iteration with the lowest
% relative residual and use that as the reported reconstruction.
mua_best = mua_est;
resid_best = inf;
iter_best = 0;

for iter = 1:maxiter
    cfg.vol = reshape(mua_est, 1, Nx, Ny, Nz);

    fprintf('\n=== iter %d ===\n', iter);
    fprintf('  forward + Jacobian (mcxlab) ...\n');
    tic;
    [detphi, ~, ~, ~, ~, Jext] = rbrunforward(cfg);
    fprintf('    done in %.2f s\n', toc);

    r = y_meas - detphi(:);
    cur_resid = norm(r) / norm(y_meas);
    resid_history(iter) = cur_resid;
    fprintf('  relative residual: %.4e\n', cur_resid);

    if (cur_resid < resid_best)
        resid_best = cur_resid;
        mua_best = mua_est;
        iter_best = iter;
    end

    if (cur_resid < 1e-4)
        fprintf('  converged.\n');
        break
    end

    % alpha scaled to the current Jacobian's response magnitude
    Jtr = rbjacop(Jext.mua);
    grad0 = Jtr(r, 'transp');
    alpha = alpha_rel * max(abs(grad0));
    fprintf('  TV-FISTA: alpha = %.3e\n', alpha);

    tic;
    [delta_mu_step, info] = rbregtv(Jext.mua, r, alpha, ...
                                    'maxit', 200, 'innerit', 10, ...
                                    'tol', 1e-6, 'verbose', 1, ...
                                    'positive', true);
    fprintf('    done in %.2f s   FISTA iters: %d   final obj: %.4e   adjoint err: %.2e\n', ...
            toc, info.iter, info.obj(end), info.adjoint_err);

    mua_est = mua_est + single(delta_mu_step);

    % keep mua physically non-negative
    mua_est = max(mua_est, single(0));
end

% adopt the best-residual iterate as the reconstruction
fprintf('\nBest iterate: iter %d with rel residual %.4e\n', iter_best, resid_best);
mua_est = mua_best;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Compare reconstruction to ground truth
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

delta_mu_rec = double(mua_est) - mua_bg;
slice_z = 15;
true_slice = delta_mu_true(:, :, slice_z);
rec_slice = delta_mu_rec(:, :, slice_z);

fprintf('\nGround-truth sum (3D): %.3e   Reconstructed sum: %.3e   ratio: %.3f\n', ...
        sum(delta_mu_true(:)), sum(delta_mu_rec(:)), ...
        sum(delta_mu_rec(:)) / max(sum(delta_mu_true(:)), eps));
fprintf('Ground-truth peak (z=%d): %.3e   Recon peak: %.3e\n', ...
        slice_z, max(true_slice(:)), max(rec_slice(:)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Plot truth vs recon and residual / objective histories
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Name', 'mcxlab iterative TV-FISTA voxel-grid mua reconstruction');
subplot(2, 2, 1);
imagesc(true_slice');
axis image;
colorbar;
set(gca, 'YDir', 'normal');
title(sprintf('ground-truth delta mua  (z=%d)', slice_z));
xlabel('x (voxel)');
ylabel('y (voxel)');

subplot(2, 2, 2);
imagesc(rec_slice');
axis image;
colorbar;
set(gca, 'YDir', 'normal');
title(sprintf('TV recon (best iter %d)  (z=%d)', iter_best, slice_z));
xlabel('x (voxel)');
ylabel('y (voxel)');

subplot(2, 2, 3);
semilogy(1:iter, resid_history(1:iter), '-o');
hold on;
semilogy(iter_best, resid_best, 'rs', 'MarkerSize', 12, 'LineWidth', 2);
hold off;
xlabel('GN iteration');
ylabel('rel residual ||y_{meas} - y_{model}|| / ||y_{meas}||');
title('residual history');
grid on;

subplot(2, 2, 4);
semilogy(info.obj, '-');
xlabel('FISTA iteration (last GN step)');
ylabel('objective  0.5||J*x-r||^2 + \alpha TV(x)');
title('FISTA convergence');
grid on;
