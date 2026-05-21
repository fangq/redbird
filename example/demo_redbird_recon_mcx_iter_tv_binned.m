%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Redbird - Iterative GN voxel-grid mua recon via mcxlab
% combining COARSE BINNING + TV-FISTA + POSITIVITY
%
% This is the "all three priors at once" version of the iterative
% recon.  It addresses the three failure modes seen in the simpler
% demos:
%
%   1. demo_redbird_recon_mcx_iter.m       (fine + LSQR):
%        single centroidal blob; LSQR's min-norm bias smears two
%        truth inclusions into one.
%
%   2. demo_redbird_recon_mcx_tv.m         (fine + TV-FISTA):
%        TV bias toward piecewise-constant solutions, but on the
%        fine grid (Nv = 1.08e5 vs Nsd = O(few hundred)) the problem
%        is still grossly underdetermined, recon amplitude lives in
%        float32 quantization noise.
%
%   3. demo_redbird_recon_mcx_iter.m       (binned + damped LSQR):
%        Tractable Nv_c, GN converges; but Tikhonov damp + min-norm
%        still gives one smeared centroidal blob with negative
%        ringing artifacts.
%
% This demo:
%   * bins the fine voxel Jacobian to a coarse reconstruction basis
%     (rbbinjac), bringing Nv_c down to ~1e4 -- ~22x smaller than
%     fine, comparable to a FEM mesh's node count.  This makes the
%     GN linearization tractable.
%   * solves each GN step with rbregtv (FISTA + isotropic 3D TV +
%     positivity).  TV biases toward piecewise-constant inclusions;
%     positivity removes the wrong-sign ambiguity that lets LSQR
%     fit MC noise.
%   * upsamples the coarse update back to the fine cfg.vol for the
%     next mcxlab forward call.
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

% ground-truth perturbation: two inclusions
delta_mu_true = zeros(Nx, Ny, Nz);
delta_mu_true(10:15, 28:32, 14:16) = 0.005;
delta_mu_true(28:32, 38:42, 14:16) = 0.008;

% sources / detectors: symmetric 5 x 5 + 5 x 5 (Nsd = 625)
[Sx, Sy] = ndgrid([10 20 30 40 50], [10 20 30 40 50]);
cfg_template.srcpos = [Sx(:), Sy(:), zeros(numel(Sx), 1)];
cfg_template.srcdir = [0  0  1];
cfg_template.srctype = 'pencil';
[Dx, Dy] = ndgrid([15 22 30 38 45], [15 22 30 38 45]);
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
cfg_truth.nphoton = 5e8;

fprintf('Synthesizing target measurements (mcxlab, true mua) ...\n');
tic;
[y_meas, ~] = rbrunforward(cfg_truth);
fprintf('  done in %.2f s   y_meas shape: %s\n', toc, mat2str(size(y_meas)));

y_meas = y_meas(:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Gauss-Newton loop with binned-TV-FISTA inner solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mua_est = single(mua_bg * ones(Nx, Ny, Nz));
cfg = cfg_template;
cfg.nphoton = 5e8;

maxiter = 5;
resid_history = zeros(maxiter, 1);

% Coarse recon basis size.  binsize = 2 keeps spatial resolution
% reasonable for resolving 5-voxel-wide truth inclusions (each spans
% ~3x3 coarse cells).  For tighter inclusions or fewer measurements,
% try binsize = 3 (more aggressive smoothing) or 1 (full fine grid).
binsize = 2;

% TV regularization weight, scaled to the current Jacobian's response
% magnitude (recomputed each GN iter).  alpha_rel ~ 1e-4 is a good
% starting point with positivity active; raise for blockier recon,
% lower for sharper edges (at the risk of letting noise back in).
alpha_rel = 1e-4;

% Stability bookkeeping: track best post-update iterate.
mua_best = mua_est;
resid_best = inf;
iter_best = 0;
info_best = struct();

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

    if (iter == 1)
        resid_best = cur_resid;
        mua_best = mua_est;
        iter_best = 0;
    elseif (cur_resid < resid_best)
        resid_best = cur_resid;
        mua_best = mua_est;
        iter_best = iter - 1;
    end

    if (cur_resid < 1e-4)
        fprintf('  converged.\n');
        break
    end

    % bin to coarse basis
    fprintf('  bin J to coarse basis (binsize = %d) ...\n', binsize);
    J_c = rbbinjac(Jext.mua, binsize);
    Nv_c = numel(J_c) / size(J_c, 4);
    fprintf('    coarse unknowns: %d   (fine: %d   ratio: %.0fx)\n', ...
            Nv_c, Nx * Ny * Nz, (Nx * Ny * Nz) / Nv_c);

    % scale alpha to the current ||J_c'r||_inf so the TV penalty
    % tracks the data-fit magnitude across GN iters
    Jc_op = rbjacop(J_c);
    grad0 = Jc_op(r, 'transp');
    alpha = alpha_rel * max(abs(grad0));

    fprintf('  TV-FISTA (alpha = %.3e, positive) ...\n', alpha);
    tic;
    [dmu_c, info] = rbregtv(J_c, r, alpha, ...
                            'maxit', 200, 'innerit', 10, ...
                            'tol', 1e-6, 'verbose', 0, ...
                            'positive', true);
    fprintf('    done in %.2f s   FISTA iters: %d   final obj: %.4e   adjoint err: %.2e\n', ...
            toc, info.iter, info.obj(end), info.adjoint_err);

    delta_mu_step = rbupsample(dmu_c, binsize, [Nx, Ny, Nz]);

    mua_est = mua_est + single(delta_mu_step);
    mua_est = max(mua_est, single(0));

    info_best = info;     % keep the last FISTA convergence trace for plotting
end

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
%%   Plot truth vs recon, residual history, FISTA trace
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Name', 'mcxlab GN voxel mua recon: binned + TV-FISTA + positivity');
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
title(sprintf('TV-binned recon (best iter %d)  (z=%d)', iter_best, slice_z));
xlabel('x (voxel)');
ylabel('y (voxel)');

subplot(2, 2, 3);
semilogy(1:iter, resid_history(1:iter), '-o');
hold on;
semilogy(max(iter_best, 1), resid_best, 'rs', 'MarkerSize', 12, 'LineWidth', 2);
hold off;
xlabel('GN iteration');
ylabel('rel residual ||y_{meas} - y_{model}|| / ||y_{meas}||');
title('residual history');
grid on;

subplot(2, 2, 4);
if (isfield(info_best, 'obj') && ~isempty(info_best.obj))
    semilogy(info_best.obj, '-');
    xlabel('FISTA iteration (last GN step)');
    ylabel('objective  0.5||J_c x - r||^2 + \alpha TV(x)');
    title('TV-FISTA convergence');
    grid on;
end
