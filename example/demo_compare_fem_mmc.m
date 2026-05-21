%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Redbird - FEM vs MMC forward and adjoint Jacobian comparison
%
% Side-by-side validation of the FEM (rbrunforward FEM path + rbjac)
% and MMC (mmclab mesh adjoint Jacobian) implementations on a single
% homogeneous slab.  RF mode (cfg.omega > 0) so we see both the real
% and imaginary parts of phi and J.
%
% Two figures are produced:
%   1. Forward fluence comparison: log10|phi| and arg(phi)
%      sliced at the y = 30 mid-plane, FEM vs MMC.
%   2. Adjoint J_mua slice for one source-detector pair: log10|J| and
%      arg(J) on the same mid-plane, FEM vs MMC.
%
% Use contourf with overlaid sparse contour lines so the iso-value
% structure is directly readable, instead of plotmesh's interpolated
% surface coloring.
%
% This file is part of Redbird URL:http://mcx.space/#redbird
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (~exist('rbrun', 'file'))
    addpath(fullfile(pwd, '../matlab'));
end
if (~exist('mmclab', 'file'))
    error('mmclab is not on the MATLAB path; this demo requires MMC');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Mesh: homogeneous 60 x 60 x 30 mm slab
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear cfg;

% Use meshgrid6 for a uniform structured tet mesh: cube/2D-slice cuts
% give regular polygons -> much smoother contour lines than the
% unstructured meshabox output.
h = 1.5;
[cfg.node, cfg.elem] = meshgrid6(0:h:60, 0:h:60, 0:h:30);
cfg.face = volface(cfg.elem);
cfg.seg  = ones(size(cfg.elem, 1), 1);

cfg.srcpos = [30 30 0];
cfg.srcdir = [0 0 1];

% Three detectors along +x.  We'll plot the Jacobian for pair 1 (S-D
% separation 10 mm) -- shallowest banana so the structure is most
% visible in a single mid-plane slice.
cfg.detpos = [20 30 0
              40 30 0
              45 30 0];
cfg.detdir = [0 0 1; 0 0 1; 0 0 1];

% Use z = boundary - epsilon so MMC's tsearchn-based detector
% localization returns a valid element index rather than NaN.
cfg.detpos(:, 3) = 0.001;
cfg.srcpos(:, 3) = 0.001;

cfg.prop  = [0 0 1 1; 0.005 1 0 1.37];
cfg.omega = 1e8 * 2 * pi;
cfg = rbmeshprep(cfg);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Run FEM forward + FEM Jacobian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cfg_fem = cfg;
for f = {'tstart', 'tend', 'tstep'}
    if (isfield(cfg_fem, f{1}))
        cfg_fem = rmfield(cfg_fem, f{1});
    end
end
[detphi_fem, phi_fem_raw] = rbrunforward(cfg_fem);
sd = rbsdmap(cfg_fem);
Jmua_fem = rbjac(sd, phi_fem_raw, cfg_fem.deldotdel, cfg_fem.elem, cfg_fem.evol, 'rfcw', 1);

% rbrunforward returns phi as a wrapped struct/map; unwrap to a plain
% (Nn x Nsrc) complex matrix
phi_fem = phi_fem_raw;
if (isstruct(phi_fem))
    phi_fem = phi_fem(1).phi;
    keys    = phi_fem.keys;
    phi_fem = phi_fem(keys{1});
end

% Jmua_fem may also be wrapped (containers.Map per wavelength)
if (isstruct(Jmua_fem))
    keys     = Jmua_fem(1).J.keys;
    Jmua_fem = Jmua_fem(1).J(keys{1});
end

fprintf('FEM phi:   %s   class=%s\n', mat2str(size(phi_fem)),   class(phi_fem));
fprintf('FEM Jmua:  %s   class=%s\n', mat2str(size(Jmua_fem)), class(Jmua_fem));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Run MMC forward + adjoint Jacobian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cfg_mmc = cfg;
cfg_mmc.nphoton  = 5e7;
cfg_mmc.tstart   = 0;
cfg_mmc.tend     = 5e-9;
cfg_mmc.tstep    = 5e-9;
cfg_mmc.gpuid    = 1;
cfg_mmc.isreflect = 1;
cfg_mmc.srcdir = [0 0 1 -inf];
% if (isfield(cfg_mmc, 'detdir'))
%     cfg_mmc = rmfield(cfg_mmc, 'detdir');
% end

[detphi_mmc, phi_mmc, ~, ~, ~, Jext_mmc] = rbrunforward(cfg_mmc);
Jmua_mmc = Jext_mmc.mua;

fprintf('MMC phi:   %s   class=%s\n', mat2str(size(phi_mmc)),  class(phi_mmc));
fprintf('MMC Jmua:  %s   class=%s\n', mat2str(size(Jmua_mmc)), class(Jmua_mmc));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Slice the tet mesh at y = 30 and interpolate to a
%%   regular grid using qmeshcut + griddata (the same
%%   pattern used in other redbird/example demos).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

phi_fem_src1 = phi_fem(:, 1);                % first (and only) source
phi_mmc_src1 = phi_mmc(:, 1);

[Xq, Zq] = meshgrid(0.5:60 - 0.5, 0.5:30 - 0.5);
slicestr = 'y=30';

% complex values: cut Re and Im parts separately, then recombine
[cutpos, cv_fem_re] = qmeshcut(cfg.elem, cfg.node, real(phi_fem_src1), slicestr);
[~,      cv_fem_im] = qmeshcut(cfg.elem, cfg.node, imag(phi_fem_src1), slicestr);
phi_fem_re_grid = griddata(cutpos(:, 1), cutpos(:, 3), cv_fem_re, Xq, Zq);
phi_fem_im_grid = griddata(cutpos(:, 1), cutpos(:, 3), cv_fem_im, Xq, Zq);
phi_fem_grid    = phi_fem_re_grid + 1i * phi_fem_im_grid;

[cutpos, cv_mmc_re] = qmeshcut(cfg.elem, cfg.node, real(phi_mmc_src1), slicestr);
[~,      cv_mmc_im] = qmeshcut(cfg.elem, cfg.node, imag(phi_mmc_src1), slicestr);
phi_mmc_re_grid = griddata(cutpos(:, 1), cutpos(:, 3), cv_mmc_re, Xq, Zq);
phi_mmc_im_grid = griddata(cutpos(:, 1), cutpos(:, 3), cv_mmc_im, Xq, Zq);
phi_mmc_grid    = phi_mmc_re_grid + 1i * phi_mmc_im_grid;

% --- Plot forward fluence comparison ---
% Use explicit contour levels so FEM (red, solid) and MMC (blue,
% dashed) lines can be overlaid on the same axes -- identical
% solutions then have lines that lie on top of each other.
mag_levels = -10:0.5:0;
arg_levels = -pi:pi / 20:pi;

figure('Name', 'FEM vs MMC forward fluence (y=30 slice, RF 70 MHz)', ...
       'Position', [60 60 900 700]);

subplot(2, 1, 1);
contour(Xq, Zq, log10(abs(phi_fem_grid) + 1e-30), mag_levels, 'r-', 'linewidth', 1.2);
hold on;
contour(Xq, Zq, log10(abs(phi_mmc_grid) + 1e-30), mag_levels, 'b--', 'linewidth', 1.2);
plot(cfg.srcpos(1, 1), cfg.srcpos(1, 3), 'kp', 'markerfacecolor', 'y', 'markersize', 14);
plot(cfg.detpos(:, 1), cfg.detpos(:, 3), 'ks', 'markerfacecolor', 'g', 'markersize', 10);
axis equal tight;
xlabel('x (mm)');
ylabel('z (mm)');
title('log_{10}|\phi|  contours: FEM (red solid) vs MMC (blue dashed),  y=30');
legend('FEM', 'MMC', 'source', 'detectors', 'location', 'southeast');

subplot(2, 1, 2);
contour(Xq, Zq, angle(phi_fem_grid), arg_levels, 'r-', 'linewidth', 1.2);
hold on;
contour(Xq, Zq, angle(phi_mmc_grid), arg_levels, 'b--', 'linewidth', 1.2);
plot(cfg.srcpos(1, 1), cfg.srcpos(1, 3), 'kp', 'markerfacecolor', 'y', 'markersize', 14);
plot(cfg.detpos(:, 1), cfg.detpos(:, 3), 'ks', 'markerfacecolor', 'g', 'markersize', 10);
axis equal tight;
xlabel('x (mm)');
ylabel('z (mm)');
title('arg(\phi)  contours: FEM (red solid) vs MMC (blue dashed),  y=30');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Adjoint J_mua for ONE source-detector pair (pair 1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% pair index k = 1 -> source 1, detector 1 (S-D separation 10 mm)
% FEM:  Jmua_fem(k, :) is the Nn-vector for pair k
% MMC:  Jmua_mmc(k, :) is the Nn-vector for pair k

pair_id = 1;

J_fem_k = Jmua_fem(pair_id, :).';
J_mmc_k = Jmua_mmc(pair_id, :).';

[cutpos, cvJ_fem_re] = qmeshcut(cfg.elem, cfg.node, real(J_fem_k), slicestr);
[~,      cvJ_fem_im] = qmeshcut(cfg.elem, cfg.node, imag(J_fem_k), slicestr);
J_fem_re_grid = griddata(cutpos(:, 1), cutpos(:, 3), cvJ_fem_re, Xq, Zq);
J_fem_im_grid = griddata(cutpos(:, 1), cutpos(:, 3), cvJ_fem_im, Xq, Zq);
J_fem_grid    = J_fem_re_grid + 1i * J_fem_im_grid;

[cutpos, cvJ_mmc_re] = qmeshcut(cfg.elem, cfg.node, real(J_mmc_k), slicestr);
[~,      cvJ_mmc_im] = qmeshcut(cfg.elem, cfg.node, imag(J_mmc_k), slicestr);
J_mmc_re_grid = griddata(cutpos(:, 1), cutpos(:, 3), cvJ_mmc_re, Xq, Zq);
J_mmc_im_grid = griddata(cutpos(:, 1), cutpos(:, 3), cvJ_mmc_im, Xq, Zq);
J_mmc_grid    = J_mmc_re_grid + 1i * J_mmc_im_grid;

% MMC's J has a known 4*Ns*Nd magnitude inflation; the SHAPE
% (banana, phase pattern) is what we compare here.  Normalize each
% so the contour LEVELS are over the same dynamic range -- this way
% same-shape distributions land on the same lines.
J_fem_scaled = J_fem_grid;
J_mmc_scaled = J_mmc_grid;

figure('Name', 'FEM vs MMC adjoint J_mua (y=30 slice, pair 1, RF 70 MHz)', ...
       'Position', [80 80 900 700]);

subplot(2, 1, 1);
contour(Xq, Zq, log10(abs(J_fem_scaled) + 1e-30), mag_levels, 'r-', 'linewidth', 1.2);
hold on;
contour(Xq, Zq, log10(abs(J_mmc_scaled) + 1e-30), mag_levels, 'b--', 'linewidth', 1.2);
plot(cfg.srcpos(1, 1), cfg.srcpos(1, 3), 'kp', 'markerfacecolor', 'y', 'markersize', 14);
plot(cfg.detpos(pair_id, 1), cfg.detpos(pair_id, 3), 'ks', 'markerfacecolor', 'g', 'markersize', 10);
axis equal tight;
xlabel('x (mm)');
ylabel('z (mm)');
title(sprintf('log_{10}|J_{\\mu_a}|  pair %d:  FEM (red) vs MMC (blue dashed)', pair_id));
legend('FEM', 'MMC', 'source', 'detector', 'location', 'southeast');

subplot(2, 1, 2);
contour(Xq, Zq, angle(J_fem_grid), arg_levels, 'r-', 'linewidth', 1.2);
hold on;
contour(Xq, Zq, angle(J_mmc_grid), arg_levels, 'b--', 'linewidth', 1.2);
plot(cfg.srcpos(1, 1), cfg.srcpos(1, 3), 'kp', 'markerfacecolor', 'y', 'markersize', 14);
plot(cfg.detpos(pair_id, 1), cfg.detpos(pair_id, 3), 'ks', 'markerfacecolor', 'g', 'markersize', 10);
axis equal tight;
xlabel('x (mm)');
ylabel('z (mm)');
title(sprintf('arg(J_{\\mu_a})  pair %d:  FEM (red) vs MMC (blue dashed)', pair_id));

% (no separate panels — overlay is what we want for direct comparison) %#ok<COMMENTNOTOK>

% --- DELETE the rest of the old 2x2 figure (replaced above) ---

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Numerical agreement summary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mask_phi = (abs(phi_fem_src1) > 1e-12) & (abs(phi_mmc_src1) > 1e-12);
ccmag = corrcoef(log10(abs(phi_fem_src1(mask_phi))), log10(abs(phi_mmc_src1(mask_phi))));
ccarg = corrcoef(angle(phi_fem_src1(mask_phi)),     angle(phi_mmc_src1(mask_phi)));
fprintf('\nForward phi agreement:\n');
fprintf('  log10|phi|  corr coef = %.4f\n', ccmag(1, 2));
fprintf('  arg(phi)    corr coef = %.4f\n', ccarg(1, 2));
fprintf('  |phi_mmc|/|phi_fem|: median=%.3f   mean=%.3f\n', ...
        median(abs(phi_mmc_src1(mask_phi)) ./ abs(phi_fem_src1(mask_phi))), ...
        mean  (abs(phi_mmc_src1(mask_phi)) ./ abs(phi_fem_src1(mask_phi))));

mask_J = (abs(J_fem_k) > 1e-12) & (abs(J_mmc_k) > 1e-12);
ccmag = corrcoef(log10(abs(J_fem_k(mask_J))), log10(abs(J_mmc_k(mask_J))));
ccarg = corrcoef(angle(J_fem_k(mask_J)),     angle(J_mmc_k(mask_J)));
fprintf('\nAdjoint J_mua agreement (pair %d):\n', pair_id);
fprintf('  log10|J|  corr coef = %.4f\n', ccmag(1, 2));
fprintf('  arg(J)    corr coef = %.4f\n', ccarg(1, 2));
fprintf('  |J_mmc|/|J_fem|: median=%.3f   mean=%.3f\n', ...
        median(abs(J_mmc_k(mask_J)) ./ abs(J_fem_k(mask_J))), ...
        mean  (abs(J_mmc_k(mask_J)) ./ abs(J_fem_k(mask_J))));
