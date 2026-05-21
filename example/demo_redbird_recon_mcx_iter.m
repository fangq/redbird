%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Redbird - Iterative Gauss-Newton voxel-grid mua recon via mcxlab
%
% Companion to demo_redbird_recon_mcx.m.  The single-step demo
% linearizes around a homogeneous background and runs ONE LSQR step;
% with the typical sparse DOT optode array, the minimum-L2-norm
% property of LSQR produces a single smooth blob even when the truth
% has multiple inclusions.
%
% This demo wraps the LSQR step in an outer Gauss-Newton loop.  At each
% iteration:
%   1. Re-run mcxlab with the CURRENT per-voxel mua estimate to obtain
%      the model prediction y_model and the adjoint Jacobian J.
%   2. Solve   J * delta_mu = (y_meas - y_model)   via rbreglsqr (LSQR).
%   3. Apply the update to cfg.vol.
%
% mcxlab consumes per-voxel mua via the MEDIA_MUA_FLOAT mode: cfg.vol
% must be a single-precision 4D array of shape (1, Nx, Ny, Nz), where
% the leading singleton axis tells mcxlab "1 property channel per
% voxel (mua)".  mus, g, n still come from cfg.prop's tissue rows.
%
% Important caveat: even with the GN outer loop, an underdetermined
% problem (Nv >> Nsd) does NOT magically resolve compact inclusions.
% The minimum-norm bias persists.  Iterating does help relative to a
% single step because:
%   - the linearization point moves with each estimate, so different J
%     "banana" shapes get sampled across iterations;
%   - the predicted model y_model from mcxlab is the true nonlinear
%     forward, so the residual r = y_meas - y_model carries information
%     the single-step linearization missed.
% For high-resolution sparse-inclusion recovery, combine this with TV /
% L1 regularization or a coarse reconstruction basis (separate work).
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

mua_bg = 0.005;          % background absorption (1/mm)
musp_bg = 1.0;            % background reduced scattering (1/mm)
n_bg = 1.37;             % refractive index

% ground-truth perturbation: two inclusions, doubled extent in each
% dimension vs the original (6x5x3 / 5x5x3) -> (12x10x6 / 10x10x6) voxels,
% giving 8x the volume per inclusion.  Bigger targets are easier for the
% min-norm LSQR to recover (more signal in r = J' * dmu_true).
delta_mu_true = zeros(Nx, Ny, Nz);
delta_mu_true(10:21, 28:37, 14:19) = 0.005;
delta_mu_true(28:37, 38:47, 14:19) = 0.008;

% sources / detectors: 5 x 5 = 25 sources + 5 x 5 = 25 detectors co-located
% on a regular 10-mm grid spanning x,y in {10, 20, 30, 40, 50} on the +z
% face (z = 0).  Matches the geometry of demo_redbird_recon_mcx.m.
% Covers x in [10, 50], y in [10, 50], which puts both truth inclusions
% (at (10..15, 28..32, 14..16) and (28..32, 38..42, 14..16)) well inside
% the banana sensitivity of multiple src-det pairs.  Nsd = 25 x 25 = 625
% measurements.
[Sx, Sy] = meshgrid(10:10:50);
[Dx, Dy] = meshgrid(10:10:50);
cfg_template.srcpos = [Sx(:), Sy(:), zeros(numel(Sx), 1)];
cfg_template.srcdir = repmat([0 0 1], numel(Sx), 1);
cfg_template.srctype = 'pencil';
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
% Fix the RNG seed so the y_meas (mcx @ mua_truth) and y0 (mcx @ mua_bg)
% forwards launch identical photon trajectories.  Their difference then
% retains only the deterministic mua-perturbation signal (this is the
% "perturbation Monte Carlo" trick: with the same seeds, photon paths
% are nearly identical, attenuation differences alone produce the
% measurable change).  Without this, the per-pair MC noise of the two
% independent runs adds in quadrature and dominates y_delta on low-
% fluence pairs.
cfg_template.seed = 12345;

Nsrc = size(cfg_template.srcpos, 1);
Ndet = size(cfg_template.detpos, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Synthesize "measurements"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Two modes (toggle via linearize_y_meas below):
%
% (A) linearize_y_meas = true  (matches demo_redbird_recon_mcx.m)
%     y_meas is built INSIDE iter 1 from that iteration's own
%     forward output as:
%         y_meas = detphi_iter1 + J_iter1' * delta_mu_true
%     so iter 1's r = J_iter1' * delta_mu_true exactly (no MC noise
%     mismatch between the linearization Jacobian and the target).
%     A single LSQR step then reproduces demo_redbird_recon_mcx.m's
%     reconstruction.  Useful only for maxiter = 1; iter 2+ would see
%     a stale y_meas built around iter 1's linearization point.
%
% (B) linearize_y_meas = false  (the original "honest" target)
%     Run mcxlab at mua_bg + delta_mu_true and use the nonlinear
%     forward output as y_meas.  This is what real measurements look
%     like.  The GN loop has nontrivial work to do across iterations
%     as the linearization point moves.

linearize_y_meas = true;

% 25 sources + 25 detectors -> 50 ExtraSrc adjoint slots; bump photon count
% so each slot still gets ~4e6 photons after mcx's uniform per-slot split.
if (linearize_y_meas)
    % y_meas built later (inside iter 1) using the iter-1 forward output.
    % Leave a placeholder so the size check below doesn't complain.
    y_meas = [];
    fprintf('Linearized y_meas mode: target measurements will be built inside iter 1.\n');
else
    mua_truth = single(mua_bg + delta_mu_true);
    cfg_truth = cfg_template;
    cfg_truth.vol = reshape(mua_truth, 1, Nx, Ny, Nz);
    cfg_truth.nphoton = 1e9;

    % CRITICAL for perturbation-MC noise cancellation: match the loop's
    % mcx launch configuration so the same-seed trick produces correlated
    % trajectories.  The loop runs with needJacobian=true -> rbrunforward
    % auto-fills cfg.detdir and uses srcid=-1, outputtype='adjoint'
    % (Ns + Nd = 50 ExtraSrc slots).  Force cfg_truth to use the SAME
    % 50-slot launch by pre-filling detdir and using srcid=-2 (Ns + Nd
    % slots, no Jacobian).  Without this match, mcx's per-thread photon
    % partition differs between the two runs and seed correlation fails.
    cfg_truth.detdir = rbgetdetdir_vol(cfg_truth);
    cfg_truth.srcid = -2;
    cfg_truth.outputtype = 'fluence';

    fprintf('Synthesizing target measurements (mcxlab nonlinear, true mua) ...\n');
    tic;
    [y_meas, ~] = rbrunforward(cfg_truth);
    fprintf('  done in %.2f s   y_meas shape: %s\n', toc, mat2str(size(y_meas)));
    y_meas = y_meas(:);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Gauss-Newton loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mua_est = single(mua_bg * ones(Nx, Ny, Nz));
cfg = cfg_template;
cfg.nphoton = 1e9;

maxiter = 3;
resid_history = zeros(maxiter, 1);

% Coarse reconstruction basis (analogue of FEM-node basis in
% mesh-based redbird).  binsize = 2 bins every 2x2x2 fine voxel block
% into one coarse-recon voxel, reducing the unknown count from
% Nx*Ny*Nz = 1.08e5 down to (Nx/2)*(Ny/2)*(Nz/2) = 1.35e4.  Combined
% with Nsd = 625 measurements, this brings the underdetermination
% ratio from ~170:1 down to ~22:1 -- the LSQR step is far less
% noise-driven and the GN linearization tracks the nonlinear forward
% much better.  For an even smaller, near-determined system try
% binsize = 4 (Nv_c ~ 1800) or binsize = 5 (Nv_c ~ 864).
%
% Damping (damp_rel below) sets the LSQR Tikhonov damp as
% damp = damp_rel * sigma_max(J_c).  Larger values shrink the LSQR
% solution toward zero; set 0 for pure min-norm.
%
% Linearized-y_meas mode disables binning and damping so the iter 1
% LSQR step reproduces demo_redbird_recon_mcx.m exactly.
if (linearize_y_meas)
    binsize = 1;
    damp_rel = 0.0;
else
    binsize = 2;
    damp_rel = 1e-2;
end

% Stability bookkeeping.  iter_best tracks the GN iterate AFTER which
% the residual was lowest (iter_best = 0 means the initial bg was
% best, i.e. every update degraded the fit -- a signal that LSQR is
% still too aggressive for the current binsize / damp).
mua_best = mua_est;
resid_best = inf;
iter_best = 0;

for iter = 1:maxiter
    % cfg.vol encoding: in linearize_y_meas mode we only run one forward
    % at the homogeneous background, so use the same single-tissue (uint8)
    % encoding as demo_redbird_recon_mcx.m -- this routes mcx through the
    % same kernel path with the same Jacobian normalization.  In iterative
    % mode the per-voxel mua_est is non-uniform, so we need MEDIA_MUA_FLOAT
    % (4D single).
    if (linearize_y_meas)
        cfg.vol = uint8(ones(Nx, Ny, Nz));
    else
        cfg.vol = reshape(mua_est, 1, Nx, Ny, Nz);
    end

    fprintf('\n=== iter %d ===\n', iter);
    fprintf('  forward + Jacobian (mcxlab) ...\n');
    tic;
    [detphi, ~, ~, ~, ~, Jext] = rbrunforward(cfg);
    fprintf('    done in %.2f s\n', toc);

    % Diagnostic: print Jacobian + residual magnitudes so we can spot any
    % normalization mismatch vs demo_redbird_recon_mcx.m.
    fprintf('    detphi  : min=%.3e  max=%.3e  ||.||=%.3e\n', ...
            min(detphi(:)), max(detphi(:)), norm(detphi(:)));
    fprintf('    Jext.mua: min=%.3e  max=%.3e  ||.||=%.3e  shape=%s\n', ...
            min(Jext.mua(:)), max(Jext.mua(:)), norm(Jext.mua(:)), ...
            mat2str(size(Jext.mua)));

    % In linearized-y_meas mode, build y_meas now using THIS iteration's
    % forward output so the MC noise cancels exactly when computing r.
    % Without this, y_meas and detphi come from independent MC runs and
    % the noise difference can swamp the J' * dmu_true signal (the cause
    % of the previously observed all-zero recovery).
    if (linearize_y_meas && iter == 1)
        J2_lin = reshape(double(Jext.mua), Nx * Ny * Nz, []);
        y_meas = detphi(:) + J2_lin.' * delta_mu_true(:);
        clear J2_lin
    end

    % Diagnostic: at iter 1 the linearization point is mua_bg (uniform),
    % so detphi == y0 (baseline forward at the homogeneous background).
    % The signed PERTURBATION caused by delta_mu_true is
    %     y_delta_mcx = y_meas - y0           (true nonlinear, from mcx)
    %     y_delta_lin = J^T * delta_mu_true   (1st-order linear prediction)
    % Plot y_delta directly (much smaller magnitude than y_meas) so the
    % linearization error is visible.  Jacobian is non-positive, so the
    % delta is negative (more mua -> less detected light).
    if (iter == 1 && ~linearize_y_meas)
        % Sanity-check J shape + magnitude vs delta_mu_true.  If
        % y_delta_lin = J' * dmu_true is much smaller than y_delta_mcx,
        % either (a) J is small in the inclusion region (mcx Jacobian
        % issue in MEDIA_MUA_FLOAT mode), (b) the reshape collapses
        % unwanted leading dimensions, or (c) delta_mu_true and J's
        % voxel-axis ordering disagree.
        fprintf('  [diag] size(Jext.mua) = %s   size(delta_mu_true) = %s\n', ...
                mat2str(size(Jext.mua)), mat2str(size(delta_mu_true)));
        J2_iter1 = reshape(double(Jext.mua), Nx * Ny * Nz, []);
        fprintf('  [diag] size(J2_iter1) = %s   Nsd_expected = %d\n', ...
                mat2str(size(J2_iter1)), Nsrc * Ndet);

        % Fraction of |J|'s L1 energy that lies inside the truth bounding
        % box, per src-det pair.  Small everywhere -> J is dim where the
        % inclusion sits.
        truth_mask = delta_mu_true ~= 0;
        truth_mask_flat = truth_mask(:);
        frac_in_truth = sum(abs(J2_iter1(truth_mask_flat, :)), 1) ./ ...
                        max(sum(abs(J2_iter1), 1), eps);
        fprintf('  [diag] truth voxels = %d   J energy fraction in truth: min=%.2e  max=%.2e  median=%.2e\n', ...
                sum(truth_mask_flat), min(frac_in_truth), max(frac_in_truth), median(frac_in_truth));

        y0 = detphi(:);                         % baseline at homogeneous bg
        y_delta_lin = J2_iter1.' * delta_mu_true(:);
        y_delta_mcx = y_meas - y0;

        % Flag the collocated (k_s == k_d) pairs in the diagnostic plot
        % so we can see how much of the spiky behavior comes from them.
        k_pair_idx = 1:Nsrc * Ndet;
        is_collocated = mod(k_pair_idx - 1, Ndet + 1) == 0;
        fprintf('  [diag] %d collocated pairs (k_s == k_d) detected\n', sum(is_collocated));

        fprintf('  [diag] ||y_delta_lin|| = %.3e   ||y_delta_mcx|| = %.3e   ratio = %.3e\n', ...
                norm(y_delta_lin), norm(y_delta_mcx), ...
                norm(y_delta_lin) / max(norm(y_delta_mcx), eps));

        % Visualize one sensitivity volume: pick the src/det closest to
        % the inclusion centroid, plot J(:,:,z=mid,pair) with the truth
        % bbox overlaid.  If the banana doesn't light up the truth region,
        % the mcx adjoint Jacobian is broken in this mode.
        [ix, iy, iz] = ind2sub(size(delta_mu_true), find(truth_mask));
        c_inc = [mean(ix), mean(iy), mean(iz)];
        [~, k_s] = min(sum((cfg_template.srcpos(:, 1:2) - c_inc(1:2)) .^ 2, 2));
        [~, k_d] = min(sum((cfg_template.detpos(:, 1:2) - c_inc(1:2)) .^ 2, 2));
        k_pair = (k_s - 1) * Ndet + k_d;        % detector-fastest order
        Jslice = squeeze(Jext.mua(:, :, round(c_inc(3)), k_pair));
        figure('Name', sprintf('J slice  pair k_s=%d -> k_d=%d', k_s, k_d));
        imagesc(Jslice');
        axis image;
        set(gca, 'YDir', 'normal');
        colorbar;
        hold on;
        rectangle('Position', [min(ix) - 0.5, min(iy) - 0.5, ...
                               max(ix) - min(ix) + 1, max(iy) - min(iy) + 1], ...
                  'EdgeColor', 'r', 'LineWidth', 2);
        plot(cfg_template.srcpos(k_s, 1), cfg_template.srcpos(k_s, 2), 'g*', 'MarkerSize', 14, 'LineWidth', 2);
        plot(cfg_template.detpos(k_d, 1), cfg_template.detpos(k_d, 2), 'wo', 'MarkerSize', 14, 'LineWidth', 2);
        hold off;
        xlabel('x (voxel)');
        ylabel('y (voxel)');
        title(sprintf('J(:,:,z=%d,pair=%d)  max|J|=%.2e\nsrc(%g,%g) -> det(%g,%g)   truth bbox in red', ...
                      round(c_inc(3)), k_pair, max(abs(Jslice(:))), ...
                      cfg_template.srcpos(k_s, 1), cfg_template.srcpos(k_s, 2), ...
                      cfg_template.detpos(k_d, 1), cfg_template.detpos(k_d, 2)));

        figure('Name', 'y_delta_lin vs y_delta_mcx @ iter 1');
        npair = numel(y_meas);

        subplot(1, 3, 1);
        plot(1:npair, y_delta_mcx, 'b.-', 'DisplayName', 'y_{\delta,mcx} = y_{mcx} - y_0');
        hold on;
        plot(1:npair, y_delta_lin, 'r.--', 'DisplayName', 'y_{\delta,lin} = J^T \delta\mu');
        hold off;
        grid on;
        xlabel('src-det pair index');
        ylabel('y_\delta (signed)');
        title('per-pair perturbation');
        legend('Location', 'best');

        subplot(1, 3, 2);
        plot(y_delta_mcx, y_delta_lin, '.');
        hold on;
        vmin = min(min(y_delta_mcx), min(y_delta_lin));
        vmax = max(max(y_delta_mcx), max(y_delta_lin));
        plot([vmin, vmax], [vmin, vmax], 'k-', 'LineWidth', 1);
        hold off;
        grid on;
        axis equal;
        xlabel('y_{\delta,mcx}');
        ylabel('y_{\delta,lin}');
        title(sprintf('scatter (||rel diff||=%.2e)', ...
                      norm(y_delta_lin - y_delta_mcx) / max(norm(y_delta_mcx), eps)));

        subplot(1, 3, 3);
        plot(1:npair, y_delta_lin - y_delta_mcx, '.-');
        grid on;
        xlabel('src-det pair index');
        ylabel('y_{\delta,lin} - y_{\delta,mcx}');
        title('linearization error per pair');

        clear J2_iter1 y0 y_delta_lin y_delta_mcx
    end

    r = y_meas - detphi(:);
    cur_resid = norm(r) / norm(y_meas);
    resid_history(iter) = cur_resid;
    fprintf('  relative residual: %.4e   ||r|| = %.3e\n', cur_resid, norm(r));

    % Track the best update that has actually been applied so far.
    % cur_resid measures mua_est going INTO this iter (the state left
    % by iter k-1's LSQR step), so the "best post-update state"
    % corresponds to iter_best = iter - 1.  iter_best = 0 means the
    % initial homogeneous mua_bg was the best (i.e., every applied
    % update degraded the fit -- a signal that LSQR is too aggressive).
    if (iter == 1)
        resid_best = cur_resid;
        mua_best = mua_est;
        iter_best = 0;
    elseif (cur_resid < resid_best)
        resid_best = cur_resid;
        mua_best = mua_est;
        iter_best = iter - 1;
    end

    % Early-exit on a small residual: ONLY check from iter >= 2, otherwise
    % we'd bail without ever applying an LSQR update.  The residual
    % entering iter 1 reflects the linearized signal magnitude (very small
    % when seed-matched perturbation MC succeeds in cancelling noise --
    % e.g. ~3e-5 ratio = ||y_delta_lin|| / ||y_meas|| -- which would falsely
    % trigger a 1e-4 absolute threshold even though the LSQR step is still
    % needed to actually recover delta_mu).  From iter 2+ the residual
    % shrinks RELATIVE to the iter-1 residual as the linearization point
    % closes in on the truth, so a stagnation check makes sense.
    if (iter > 1 && cur_resid < 1e-4)
        fprintf('  converged.\n');
        break
    end

    % Bin the fine voxel Jacobian to the coarse recon basis, solve the
    % SMALL (~Nv_c) GN step there with damped LSQR, then upsample the
    % coarse update back to the fine grid for the next forward call.
    fprintf('  bin J to coarse basis (binsize = %d) ...\n', binsize);
    J_c = rbbinjac(Jext.mua, binsize);
    Nv_c = numel(J_c) / size(J_c, 4);
    fprintf('    coarse unknowns: %d   (fine: %d   ratio: %.0fx)\n', ...
            Nv_c, Nx * Ny * Nz, (Nx * Ny * Nz) / Nv_c);

    % Adaptive Tikhonov damp scaled to sigma_max(J_c).  Canonical
    % choice: damp = damp_rel * sigma_max regularizes singular values
    % below damp_rel * sigma_max and lets the rest pass through.
    %   damp_rel ~ 1e-3  : light damp, keeps most spectrum
    %   damp_rel ~ 1e-2  : balanced for typical DOT scales
    %   damp_rel ~ 1e-1  : heavy damp, smoother recon, very stable
    %
    % sigma_max comes from 10 power-iter steps on J_c^T J_c, cheap
    % since Nv_c << Nv_fine.
    afun_c = rbjacop(J_c);
    v_pow = randn(Nv_c, 1);
    v_pow = v_pow / norm(v_pow);
    for kpow = 1:10
        w_pow = afun_c(v_pow, 'notransp');
        v_pow = afun_c(w_pow, 'transp');
        nrm_pow = norm(v_pow);
        if (nrm_pow < eps)
            break
        end
        v_pow = v_pow / nrm_pow;
    end
    sigma_max = norm(afun_c(v_pow, 'notransp'));

    % damp_rel is set once before the loop (see binsize/damp_rel block);
    % keep it stable across GN iterations so the regularization strength
    % doesn't drift with sigma_max.
    damp_now = damp_rel * sigma_max;

    fprintf('  sigma_max(J_c) = %.3e   damp = %.3e (damp_rel = %.0e)\n', ...
            sigma_max, damp_now, damp_rel);
    fprintf('  LSQR ...\n');
    tic;
    [dmu_c, info] = rbreglsqr(J_c, r, ...
                              'maxit', 100, 'tol', 1e-8, 'damp', damp_now);
    fprintf('    done in %.2f s   LSQR iters: %d   relres: %.3e   adjoint err: %.2e\n', ...
            toc, info.iter, info.relres, info.adjoint_err);
    fprintf('    dmu_c  : ||.||=%.3e   max=%.3e   min=%.3e   shape=%s\n', ...
            norm(dmu_c(:)), max(dmu_c(:)), min(dmu_c(:)), mat2str(size(dmu_c)));

    delta_mu_step = rbupsample(dmu_c, binsize, [Nx, Ny, Nz]);
    fprintf('    delta_mu_step (fine): ||.||=%.3e   max=%.3e   min=%.3e\n', ...
            norm(delta_mu_step(:)), max(delta_mu_step(:)), min(delta_mu_step(:)));

    mua_est_pre_clip = mua_est + single(delta_mu_step);
    mua_est = max(mua_est_pre_clip, single(0));
    n_clipped = sum(mua_est_pre_clip(:) < 0);
    fprintf('    mua_est after update: ||.||=%.3e   range=[%.3e, %.3e]   clipped to 0: %d voxels\n', ...
            norm(double(mua_est(:))), min(mua_est(:)), max(mua_est(:)), n_clipped);
end

% Report the best in-loop residual (always measured BEFORE the iter k LSQR
% update is applied; iter_best = 0 means the initial homogeneous bg gave
% the smallest residual among the iterates examined).  Do NOT overwrite
% mua_est: the plot below should show the actual final iterate, otherwise
% maxiter = 1 always displays the homogeneous initial state.
fprintf('\nBest in-loop residual: %.4e at iter %d\n', resid_best, iter_best);

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
%%   Plot truth vs recon and residual history
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Name', 'mcxlab iterative GN voxel-grid mua reconstruction');
subplot(1, 3, 1);
imagesc(true_slice');
axis image;
colorbar;
set(gca, 'YDir', 'normal');
title(sprintf('ground-truth delta mua  (z=%d)', slice_z));
xlabel('x (voxel)');
ylabel('y (voxel)');

subplot(1, 3, 2);
imagesc(rec_slice');
axis image;
colorbar;
set(gca, 'YDir', 'normal');
title(sprintf('recovered (final iter %d)  (z=%d)', iter, slice_z));
xlabel('x (voxel)');
ylabel('y (voxel)');

subplot(1, 3, 3);
semilogy(1:iter, resid_history(1:iter), '-o');
hold on;
semilogy(iter_best, resid_best, 'rs', 'MarkerSize', 12, 'LineWidth', 2);
hold off;
xlabel('GN iteration');
ylabel('relative residual ||y_{meas} - y_{model}|| / ||y_{meas}||');
title('residual history');
grid on;
