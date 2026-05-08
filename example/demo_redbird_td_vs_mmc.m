%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Redbird time-domain forward solver vs. mmclab Monte Carlo
%
% Both solvers run on the *same tetrahedral mesh* with the *same
% source/detector geometry*, the *same optical properties*, and the *same
% time grid* (cfg.tstart/cfg.tstep/cfg.tend). The script compares the
% temporal point-spread function (TPSF) at each detector and reports a
% normalized agreement metric.
%
% Redbird uses an implicit Crank-Nicolson FEM time-stepper. mmclab uses
% mesh-based Monte Carlo with the requested time gates.
%
% This file is part of Redbird URL:http://mcx.space/#redbird
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(fullfile(pwd, '../matlab'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% shared mesh and optical properties
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear cfg mcfg;

[cfg.node, cfg.face, cfg.elem] = meshabox([0 0 0], [60 60 30], 3);
cfg.seg = ones(size(cfg.elem, 1), 1);

% optical properties: [mua mus g n] in 1/mm, 1/mm, scalar, scalar
% (typical breast/brain tissue parameters at 760 nm)
cfg.prop = [0 0 1 1
            0.01 1.0 0 1.37];

cfg.srcpos = [30 30 0];
cfg.srcdir = [0 0 1];

% three detectors at increasing source-detector separation
cfg.detpos = [30 30 30
              40 30 30
              50 30 30];
cfg.detdir = [0 0 -1];

% time grid (seconds)
cfg.tstart = 0;
cfg.tstep  = 50e-12;        % 50 ps
cfg.tend   = 4e-9;          % 4 ns

cfg.omega = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Redbird time-domain forward (Crank-Nicolson)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cfg = rbmeshprep(cfg);

fprintf(1, 'Running Redbird time-domain forward solver...\n');
tic;
[detphi_rb, ~] = rbrunforward(cfg);
toc;

% detphi_rb has shape (Ndet x Nsrc x Nt)
tvec = cfg.tstart:cfg.tstep:cfg.tend;
nt = length(tvec);
ndet = size(cfg.detpos, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% mmclab time-domain Monte Carlo on the same mesh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (~exist('mmclab', 'file'))
    error('mmclab is not on the MATLAB path; this comparison requires MMC');
end

mcfg.nphoton = 1e7;             % bump up for less Monte Carlo noise
mcfg.node = cfg.node;
mcfg.elem = cfg.elem;
mcfg.elemprop = cfg.seg;
mcfg.prop = cfg.prop;
mcfg.srcpos = cfg.srcpos;
mcfg.srcdir = cfg.srcdir;
mcfg.tstart = cfg.tstart;
mcfg.tstep  = cfg.tstep;
mcfg.tend   = cfg.tend;
mcfg.method = 'elem';           % branchless Badouel (best for TD)
mcfg.isreflect = 1;
mcfg.basisorder = 1;            % node-based output (matches FEM linear basis)
mcfg.seed = 19790521;
mcfg.outputtype = 'fluence';

% mmclab does not use cfg.detpos to build a TPSF directly; we'll interpolate
% the volumetric fluence at the detector positions afterwards.
fprintf(1, 'Running mmclab Monte Carlo time-domain...\n');
tic;
flux = mmclab(mcfg);
toc;

% flux.data has shape (Nnode x Nt). MC time-gate i covers
% [tstart+(i-1)*tstep, tstart+i*tstep], so the centred gate-time is
% tstart+(i-0.5)*tstep. For comparison with Redbird's instantaneous detphi
% sampled at the gate left-edge tvec(i+1) for i>=1, we evaluate MMC at the
% gate centres and compare in shape (peak-normalized).
mc_nt = size(flux.data, 2);
mc_tvec = mcfg.tstart + ((1:mc_nt) - 0.5) * mcfg.tstep;

% interpolate MMC fluence at each detector position via tsearchn + bary
[loc_det, bary_det] = tsearchn(cfg.node(:, 1:3), cfg.elem(:, 1:4), cfg.detpos(:, 1:3));
detphi_mc = zeros(ndet, mc_nt);
for i = 1:ndet
    if (~isnan(loc_det(i)))
        nodes_i = cfg.elem(loc_det(i), 1:4);
        for k = 1:mc_nt
            detphi_mc(i, k) = bary_det(i, :) * flux.data(nodes_i, k);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% align time bases & compute agreement metrics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Redbird detphi indexed by tvec(1)..tvec(end). MC detphi indexed by gate
% centres mc_tvec. To compare, sample Redbird at the same gate centres.
detphi_rb_at_gates = zeros(ndet, mc_nt);
for i = 1:ndet
    detphi_rb_at_gates(i, :) = interp1(tvec(:), squeeze(detphi_rb(i, 1, :)), mc_tvec(:), 'linear', 0);
end

% peak-normalize each curve so the comparison is unaffected by absolute
% normalization conventions (MMC: per launched photon; Redbird: per unit
% FEM source amplitude). Shape agreement is what matters.
fprintf(1, '\n%-12s  %-12s  %-12s  %-12s\n', 'detector', 'sd-sep (mm)', 'corr coef', 'L2 err (norm)');
fprintf(1, '%s\n', repmat('-', 1, 60));
for i = 1:ndet
    rb_i = detphi_rb_at_gates(i, :);
    mc_i = detphi_mc(i, :);
    rb_norm = rb_i / max(abs(rb_i));
    mc_norm = mc_i / max(abs(mc_i));
    cc = corrcoef(rb_norm(:), mc_norm(:));
    cc = cc(1, 2);
    l2err = norm(rb_norm - mc_norm) / norm(mc_norm);
    sdsep = norm(cfg.detpos(i, 1:3) - cfg.srcpos);
    fprintf(1, '%-12d  %-12.2f  %-12.4f  %-12.4f\n', i, sdsep, cc, l2err);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot TPSFs side-by-side (peak-normalized for shape comparison)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Name', 'Redbird TD vs mmclab TD - TPSF comparison', 'Position', [100 100 1100 360]);
for i = 1:ndet
    subplot(1, ndet, i);
    rb_i = squeeze(detphi_rb(i, 1, :));
    mc_i = detphi_mc(i, :);

    % peak-normalize for shape comparison
    rb_n = rb_i / max(abs(rb_i));
    mc_n = mc_i / max(abs(mc_i));

    semilogy(tvec * 1e9, abs(rb_n), 'r-', 'LineWidth', 2);
    hold on;
    semilogy(mc_tvec * 1e9, abs(mc_n), 'b--', 'LineWidth', 1.5);
    grid on;
    xlabel('time (ns)');
    ylabel('|\Phi(t)| / max');
    sdsep = norm(cfg.detpos(i, 1:3) - cfg.srcpos);
    title(sprintf('det %d, SD=%.1f mm', i, sdsep));
    legend({'Redbird (CN)', 'mmclab (MC)'}, 'Location', 'southwest');
    ylim([1e-4 2]);
    xlim([0 cfg.tend * 1e9]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot temporal centre-of-mass (mean-time) per detector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% the mean photon arrival time <t> = int t*phi(t) dt / int phi(t) dt is a
% standard, normalization-free observable that should agree closely between
% the two solvers if the diffusion physics is captured correctly.
mean_t_rb = zeros(ndet, 1);
mean_t_mc = zeros(ndet, 1);
for i = 1:ndet
    rb_i = max(0, squeeze(detphi_rb(i, 1, :)));
    mc_i = max(0, detphi_mc(i, :)');
    mean_t_rb(i) = sum(tvec(:) .* rb_i) / sum(rb_i);
    mean_t_mc(i) = sum(mc_tvec(:) .* mc_i) / sum(mc_i);
end

fprintf(1, '\nmean photon arrival time <t> per detector:\n');
fprintf(1, '%-12s  %-14s  %-14s  %-12s\n', 'detector', 'redbird (ns)', 'mmclab (ns)', 'rel err');
fprintf(1, '%s\n', repmat('-', 1, 60));
for i = 1:ndet
    relerr = abs(mean_t_rb(i) - mean_t_mc(i)) / mean_t_mc(i);
    fprintf(1, '%-12d  %-14.4f  %-14.4f  %-12.4f\n', ...
            i, mean_t_rb(i) * 1e9, mean_t_mc(i) * 1e9, relerr);
end
