% TEST_ISRATIO  Validate the new recon.isratio flag in rbrunrecon on a small mesh.
% Claim: reconstructing from the ratio I/I0 with isratio=1 reproduces reconstructing
% from the absolute measurement I, when the baseline model equals I0 (verified for
% Born reform='real' and Rytov reform='logphase').
%
% NOTE: cfg.prop / recon.prop / the data are containers.Map HANDLES, and rbrun (and
% the isratio patch) mutate them in place, so every rbrun call is given a FRESHLY
% built cfg/recon and a deep-copied data Map.
% prepend THIS repo's matlab/ unconditionally so it shadows any other redbird copy
% that may already be on the MATLAB path (e.g. a sibling 'redbird' checkout)
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'matlab'));
if (~exist('meshabox', 'file'))
    error('iso2mesh is required on the MATLAB path to run this test');
end

% --- perturbed "absolute" measurement (sphere inclusion) ---
clear cfg0;
s0 = [70, 50, 20];
rs = 8;
[nobbx, fcbbx] = meshabox([40 0 0], [160, 120, 60], 10);
[nosp, fcsp] = meshasphere(s0, rs, 1);
[no, fc] = mergemesh(nobbx, fcbbx, nosp, fcsp);
[cfg0.node, cfg0.elem] = s2m(no, fc(:, 1:3), 1, 40, 'tetgen', [41 1 1; s0]);
cfg0.seg = cfg0.elem(:, 5);
[xi, yi] = meshgrid(60:20:140, 20:20:100);
cfg0.srcpos = [xi(:), yi(:), zeros(numel(yi), 1)];
cfg0.srcdir = [0 0 1];
cfg0.detpos = [xi(:), yi(:), 60 * ones(numel(yi), 1)];
cfg0.detdir = [0 0 -1];
cfg0.param = struct('hbo', [10 20], 'hbr', [4 8], 'scatamp500', [1 1], 'scatpow500', [1.2 1.2]);
cfg0.prop = containers.Map({'690', '830'}, ...
                           {[0 0 1 1; 0 1 0 1.37; 0 1 0 1.37], [0 0 1 1; 0 0.8 0 1.37; 0 0.8 0 1.37]});
cfg0.omega = 0;
cfg0 = rbmeshprep(cfg0);
detphi_pert = rbrun(cfg0);

% --- fixed forward + recon geometry ---
G.srcpos = cfg0.srcpos;
G.srcdir = cfg0.srcdir;
G.detpos = cfg0.detpos;
G.detdir = cfg0.detdir;
[G.node, ~, G.elem] = meshabox([40 0 0], [160, 120, 60], 10);
[R.node, ~, R.elem] = meshabox([40 0 0], [160, 120, 60], 25);
R.elem = R.elem(:, 1:4);
[R.mapid, R.mapweight] = tsearchn(R.node, R.elem, G.node);
sd = rbsdmap(buildcfg(G));

for reform = {'real', 'logphase'}
    rf = reform{1};
    detphi_base = rbrun(buildcfg(G));                       % I0 (= model at iter 1)
    ratio = containers.Map();
    for w = detphi_pert.keys
        ratio(w{1}) = detphi_pert(w{1}) ./ detphi_base(w{1});
    end

    abs_rec = rbrun(buildcfg(G), buildrecon(R), copymap(detphi_pert), sd, ...
                    'mode', 'image', 'maxiter', 1, 'lambda', 1e-2, 'reform', rf);
    reconr = buildrecon(R);
    reconr.isratio = 1;                 % the new flag is a recon field
    rat_rec = rbrun(buildcfg(G), reconr, ratio, sd, ...
                    'mode', 'image', 'maxiter', 1, 'lambda', 1e-2, 'reform', rf);
    da = abs_rec.param.hbo(:);
    dr = rat_rec.param.hbo(:);
    fprintf('reform=%-9s  max|abs-ratio| HbO = %.3e   rel = %.3e   corr = %.6f\n', ...
            rf, max(abs(da - dr)), max(abs(da - dr)) / max(abs(da - 10)), corr(da, dr));
end
exit;

% ---- helpers (fresh handles each call) ----
function cfg = buildcfg(G)
    cfg.node = G.node;
    cfg.elem = [G.elem(:, 1:4), ones(size(G.elem, 1), 1)];
    cfg.seg = cfg.elem(:, 5);
    cfg.srcpos = G.srcpos;
    cfg.srcdir = G.srcdir;
    cfg.detpos = G.detpos;
    cfg.detdir = G.detdir;
    cfg.param = struct('hbo', 10, 'hbr', 4);                    % fixed scattering via prop
    cfg.prop = containers.Map({'690', '830'}, ...
                              {[0 0 1 1; 0 1 0 1.37], [0 0 1 1; 0 0.8 0 1.37]});
    cfg.omega = 0;
    cfg = rbmeshprep(cfg);
end

function recon = buildrecon(R)
    recon.node = R.node;
    recon.elem = R.elem;
    recon.mapid = R.mapid;
    recon.mapweight = R.mapweight;
    recon.bulk = struct('hbo', 10, 'hbr', 4);
    recon.param = struct('hbo', 10, 'hbr', 4);
    recon.prop = containers.Map({'690', '830'}, {[], []});
end

function m2 = copymap(m)
    m2 = containers.Map(m.keys, m.values);
end
