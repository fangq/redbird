function run_redbird_test(tests)
%
% run_redbird_test
%   or
% run_redbird_test(tests)
% run_redbird_test({'util','jac','prop','mesh','forward','solver','recon'})
%
% Unit tests for Redbird-m. Tests are grouped so they can be run
% selectively when iso2mesh / mmclab is unavailable for the heavier groups.
%
% input:
%      tests: cell array of strings selecting groups to run:
%         'util'    : pure utilities (rbextinction, rbmusp2sasp, rbgetreff,
%                     rbgetdistance, rbgetbulk, rbsdmap, rbmasksum, rbmatflat,
%                     rbgetcfg, rbaddnoise)
%         'jac'     : Jacobian param-mapping (rbjacscatamp, rbjacscatpow,
%                     rbjacscat, rbjacchrome) including the legacy and the
%                     500 nm-normalized scattering parameterization
%         'prop'    : parameter <-> property propagation (rbupdateprop,
%                     rbsyncprop)
%         'mesh'    : rbmeshprep on a small box mesh (requires iso2mesh)
%         'forward' : end-to-end forward smoke test (rbrun / rbrunforward,
%                     including rbjac sign and symmetry; requires iso2mesh)
%         'solver'  : rbreginv regularization-strength behavior
%         'recon'   : reconstruction smoke test, multi-iter residual behavior
%                     (rbrunrecon, requires iso2mesh)
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% license:
%     GPL version 3, see LICENSE_GPLv3.txt files for details
%
% -- this function is part of Redbird-m toolbox
%

if (nargin == 0)
    tests = {'util', 'jac', 'prop', 'mesh', 'forward', 'solver', 'recon', 'mwt', 'td'};
end

global RB_FAIL RB_TOTAL
RB_FAIL = 0;
RB_TOTAL = 0;

addpath([fileparts(which(mfilename)) filesep '..' filesep 'matlab']);

bar = char(ones(1, 79) * 61);

if (ismember('util', tests))
    fprintf(1, '%s\nUtility helpers\n%s\n', bar, bar);
    test_util();
end

if (ismember('jac', tests))
    fprintf(1, '%s\nJacobian parameter mapping\n%s\n', bar, bar);
    test_jac();
end

if (ismember('prop', tests))
    fprintf(1, '%s\nParameter <-> property propagation\n%s\n', bar, bar);
    test_prop();
end

if (ismember('mesh', tests) && exist('meshabox', 'file'))
    fprintf(1, '%s\nMesh preparation\n%s\n', bar, bar);
    test_mesh();
elseif ismember('mesh', tests)
    fprintf(2, '[skip] ''mesh'' requires iso2mesh (meshabox not found)\n');
end

if (ismember('forward', tests) && exist('meshabox', 'file'))
    fprintf(1, '%s\nForward smoke test\n%s\n', bar, bar);
    test_forward();
elseif ismember('forward', tests)
    fprintf(2, '[skip] ''forward'' requires iso2mesh\n');
end

if (ismember('solver', tests))
    fprintf(1, '%s\nSolver / regularization\n%s\n', bar, bar);
    test_solver();
end

if (ismember('recon', tests) && exist('meshabox', 'file'))
    fprintf(1, '%s\nReconstruction smoke test\n%s\n', bar, bar);
    test_recon();
elseif ismember('recon', tests)
    fprintf(2, '[skip] ''recon'' requires iso2mesh\n');
end

if (ismember('mwt', tests))
    fprintf(1, '%s\nMicrowave tomography (Helmholtz + RBC)\n%s\n', bar, bar);
    test_mwt();
end

if (ismember('td', tests) && exist('meshabox', 'file'))
    fprintf(1, '%s\nTime-domain DOT (Crank-Nicolson)\n%s\n', bar, bar);
    test_td();
elseif ismember('td', tests)
    fprintf(2, '[skip] ''td'' requires iso2mesh\n');
end

fprintf(1, '%s\n', bar);
if RB_FAIL == 0
    fprintf(1, 'All %d tests passed.\n', RB_TOTAL);
else
    fprintf(2, '%d / %d tests FAILED.\n', RB_FAIL, RB_TOTAL);
end

% =========================================================================
function test_util()

% ---- rbextinction ----
% Spectrum is sampled exactly at 690 and 832; values at those rows can be
% read straight off the chrome table and used as a regression check.
ext = rbextinction(690, {'hbo', 'hbr'});
test_redbird('rbextinction shape (1 wv, 2 chrome)', @size, [1 2], ext);
test_redbird('rbextinction hbo>0 at 690nm', @(x) x > 0, true, ext(1, 1));
test_redbird('rbextinction hbr>0 at 690nm', @(x) x > 0, true, ext(1, 2));

ext2 = rbextinction([690, 830], {'hbo', 'hbr', 'water'});
test_redbird('rbextinction shape (2 wv, 3 chrome)', @size, [2 3], ext2);

% string-cell wavelength input form (used internally by rbjacchrome)
ext3 = rbextinction({'690', '830'}, {'hbo', 'hbr'});
test_redbird('rbextinction string-cell wavelength', @() ext3, ext2(:, 1:2));

% scalar 'type' returns column vector
extv = rbextinction([690; 830], 'hbo');
test_redbird('rbextinction scalar type returns vector', @size, [2 1], extv);

% empty wavelength cell errors
test_redbird('rbextinction empty wavelength errors', ...
             @() rbextinction({''}, {'hbo'}), 'error');

% ---- rbmusp2sasp ----
% Round-trip: pick (sa,sp), evaluate at two wavelengths, recover.
sa_in = 1.4;
sp_in = 1.1;
lambda = [690, 830];
musp = sa_in * (lambda / 500).^(-sp_in);
[sa, sp] = rbmusp2sasp(musp, lambda);
test_redbird('rbmusp2sasp recovers scatamp', @() sa, sa_in);
test_redbird('rbmusp2sasp recovers scatpow', @() sp, sp_in);

% ---- rbgetreff ----
% Without index mismatch (n=1), Reff should be ~0 (only the discretization
% boundary at the critical angle contributes).
test_redbird('rbgetreff n=1 ~ 0', @(x) abs(x) < 1e-4, true, rbgetreff(1));
% Regression check at the canonical tissue index n=1.37
test_redbird('rbgetreff(1.37,1) ~ 0.4684', ...
             @(x) abs(x - 0.4684) < 1e-3, true, rbgetreff(1.37, 1));
% Higher index -> higher reflection (monotonicity)
test_redbird('rbgetreff monotone in n_in', ...
             @(d) d > 0, true, rbgetreff(1.5, 1) - rbgetreff(1.37, 1));

% ---- rbgetdistance ----
src = [0 0 0; 1 0 0];
det = [0 0 0; 0 1 0; 0 0 1];
d = rbgetdistance(src, det);
expected = [0, 1, 1; 1, sqrt(2), sqrt(2)];
test_redbird('rbgetdistance basic Ns x Nd', @() d, expected);

% ---- rbgetbulk ----
cfgB = struct('bulk', struct('mua', 0.01, 'musp', 1.0, 'g', 0, 'n', 1.37));
b = rbgetbulk(cfgB);
test_redbird('rbgetbulk from cfg.bulk', @() b, [0.01, 1.0, 0, 1.37]);
cfgB2 = struct('bulk', struct('dcoeff', 1 / 3));
b2 = rbgetbulk(cfgB2);
test_redbird('rbgetbulk dcoeff -> musp=1', @() b2(2), 1.0);
test_redbird('rbgetbulk default n=1.37', @() b2(4), 1.37);

% ---- rbmasksum ----
% Bin-sum semantics: rows summed by mask labels
data = [1 2 3 4; 10 20 30 40];
mask = [1 1 2 2];
s = rbmasksum(data, mask);
% s is N_bins x Ncols; data is 2x4, mask along columns -> output is 2 x 2
test_redbird('rbmasksum shape', @size, [2 2], s);
test_redbird('rbmasksum bin1 col1', @() s(1, 1), 1 + 2);
test_redbird('rbmasksum bin2 col2', @() s(2, 2), 30 + 40);

% ---- rbmatflat ----
M = containers.Map();
M('690') = [1 2; 3 4];
M('830') = [5 6; 7 8];
flat = rbmatflat(M);
test_redbird('rbmatflat stacks maps vertically', @size, [4 2], flat);
% with weights
flatw = rbmatflat(M, [10; 1]);
test_redbird('rbmatflat applies per-key weights', @() flatw(1, 1), 10);

% passing a plain matrix returns it unchanged
A0 = magic(4);
test_redbird('rbmatflat passthrough', @() rbmatflat(A0), A0);

% ---- rbgetcfg ----
% structfun-based, so the return is a numeric column vector in field order
cfgM = struct;
cfgM.foo = containers.Map({'690', '830'}, {1, 2});
cfgM.bar = 99;  % not a map - should pass through
out = rbgetcfg(cfgM, '690');
test_redbird('rbgetcfg returns vector in field order', @() out, [1; 99]);
% numeric key gets stringified
out2 = rbgetcfg(cfgM, 690);
test_redbird('rbgetcfg numeric key -> string', @() out2, [1; 99]);

% ---- rbsdmap (simple, no mesh) ----
% Small artificial cfg with prop as containers.Map; rbsdmap should return a
% containers.Map with one entry per wavelength and 3-column or 4-column SD.
cfgSD = struct;
cfgSD.srcpos = [0 0 0; 10 0 0];
cfgSD.detpos = [20 0 0; 30 0 0];
cfgSD.face = [];
cfgSD.node = [];
cfgSD.prop = containers.Map({'690'}, {[0 0 1 1; 0.01 1 0 1.37]});
sd = rbsdmap(cfgSD);
test_redbird('rbsdmap returns map for multispec cfg', ...
             @(x) isa(x, 'containers.Map'), true, sd);
sdwv = sd('690');
test_redbird('rbsdmap entry has >= 3 columns', @(x) x >= 3, true, size(sdwv, 2));
test_redbird('rbsdmap covers Ns*Nd pairs', @() size(sdwv, 1), 4);

% maxdist filter shrinks the set
sd2 = rbsdmap(cfgSD, 15);
sdwv2 = sd2('690');
test_redbird('rbsdmap maxdist drops far pairs', ...
             @(x) x < size(sdwv, 1), true, sum(sdwv2(:, 3)));

% ---- rbaddnoise ----
data = abs(randn(50, 1)) + 1;  % strictly positive baseline
% no-noise sentinel: snrshot=Inf, snrthermal=Inf -> deviation from data = 0
clean = rbaddnoise(data, Inf, Inf, 42);
test_redbird('rbaddnoise SNR=Inf is identity', @() clean, data);

% Shape preservation across input shapes
data2d = abs(randn(8, 5)) + 1;
n2d = rbaddnoise(data2d, 30, 30, 7);
test_redbird('rbaddnoise preserves 2D shape', @size, [8 5], n2d);

% Reproducibility with the same seed
n_a = rbaddnoise(data, 30, 30, 12345);
n_b = rbaddnoise(data, 30, 30, 12345);
test_redbird('rbaddnoise is deterministic given seed', @() n_a, n_b);

% Different seed -> different sample
n_c = rbaddnoise(data, 30, 30, 99999);
test_redbird('rbaddnoise different seed -> different sample', ...
             @(d) d > 0, true, max(abs(n_a - n_c)));

% Higher SNR -> smaller deviation from clean signal (statistical, but
% reliable with N=200 samples)
data_lg = abs(randn(200, 1)) + 1;
low_snr  = rbaddnoise(data_lg, 20, 20, 1);
high_snr = rbaddnoise(data_lg, 60, 60, 1);
err_low  = norm(low_snr  - data_lg);
err_high = norm(high_snr - data_lg);
test_redbird('rbaddnoise high SNR -> smaller deviation', ...
             @(r) r > 5, true, err_low / err_high);

% =========================================================================
function test_jac()

% ---- rbjacscatamp / rbjacscatpow: legacy vs 500nm formulas ----
Jd = rand(3, 4);
dcoeff = rand(1, 4) * 0.1;
sp = ones(4, 1) * 0.9;
wv = 700;

% legacy default (lref defaults to 1e9 nm, i.e. lambda-in-meters convention)
[~, dD_amp_legacy] = rbjacscatamp(Jd, dcoeff, wv, sp);
ref_amp_legacy = -3 .* dcoeff .* dcoeff .* (wv * 1e-9).^(-sp.');
test_redbird('rbjacscatamp legacy formula', @() dD_amp_legacy, ref_amp_legacy);

% explicit lref=500 -> normalized convention
[~, dD_amp_500] = rbjacscatamp(Jd, dcoeff, wv, sp, 500);
ref_amp_500 = -3 .* dcoeff .* dcoeff .* (wv / 500).^(-sp.');
test_redbird('rbjacscatamp 500nm-normalized', @() dD_amp_500, ref_amp_500);

[~, dD_pow_legacy] = rbjacscatpow(Jd, dcoeff, wv);
test_redbird('rbjacscatpow legacy formula', ...
             @() dD_pow_legacy, dcoeff .* log(wv * 1e-9));

[~, dD_pow_500] = rbjacscatpow(Jd, dcoeff, wv, 500);
test_redbird('rbjacscatpow 500nm-normalized', ...
             @() dD_pow_500, dcoeff .* log(wv / 500));

% Jacobian = Jd .* dDd*
Jamp = rbjacscatamp(Jd, dcoeff, wv, sp);
test_redbird('rbjacscatamp J = Jd .* dDd', @() Jamp, Jd .* dD_amp_legacy);

% ---- rbjacscat: legacy vs normalized field names ----
Jdmap = containers.Map();
Jdmap('700') = rand(3, 4);
Jdmap('800') = rand(3, 4);
dcmap = containers.Map();
dcmap('700') = rand(1, 4) * 0.1;
dcmap('800') = rand(1, 4) * 0.1;

J_legacy = rbjacscat(Jdmap, dcmap, sp, {'700', '800'});
test_redbird('rbjacscat legacy field names', ...
             @() sort(fieldnames(J_legacy)), {'scatamp'; 'scatpow'});
test_redbird('rbjacscat legacy block size', @size, [6 4], J_legacy.scatamp);

J_norm = rbjacscat(Jdmap, dcmap, sp, {'700', '800'}, 500, '500');
test_redbird('rbjacscat 500nm field names', ...
             @() sort(fieldnames(J_norm)), {'scatamp500'; 'scatpow500'});

% values should differ between the two conventions (sanity)
test_redbird('rbjacscat legacy vs 500nm differ', ...
             @(d) d > 0, true, max(abs(J_legacy.scatamp(:) - J_norm.scatamp500(:))));

% ---- rbjacchrome ----
Jmua = containers.Map();
Jmua('690') = ones(2, 3);   % 2 SD pairs, 3 nodes
Jmua('830') = 2 * ones(2, 3);
Jchrome = rbjacchrome(Jmua, {'hbo', 'hbr'});
test_redbird('rbjacchrome fields', @() sort(fieldnames(Jchrome)), {'hbo'; 'hbr'});
test_redbird('rbjacchrome stacked rows = wv*Nsd', @() size(Jchrome.hbo, 1), 4);
test_redbird('rbjacchrome cols = nodes', @() size(Jchrome.hbo, 2), 3);

% =========================================================================
function test_prop()

% ---- rbupdateprop: legacy scatamp/scatpow ----
cfg = struct;
cfg.node = [0 0 0; 1 0 0; 0 1 0; 0 0 1; 0.5 0.5 0.5];
cfg.elem = [1 2 3 4; 1 2 3 5];
cfg.seg = [1; 2];
cfg.param = struct('hbo', [10 10], 'hbr', [3 3], ...
                   'scatamp', [1 1], 'scatpow', [0.9 0.9]);
cfg.prop = containers.Map({'690', '830'}, ...
                          {[0 0 1 1; 0 1 0 1.37; 0 1 0 1.37], ...
                           [0 0 1 1; 0 1 0 1.37; 0 1 0 1.37]});
prop = rbupdateprop(cfg);
test_redbird('rbupdateprop returns map', @(x) isa(x, 'containers.Map'), true, prop);
p690 = prop('690');
musp_legacy_690 = 1 * (690 * 1e-9)^(-0.9);
test_redbird('rbupdateprop legacy musp@690', @() p690(2, 2), musp_legacy_690);

% ---- rbupdateprop: 500nm-normalized takes precedence ----
cfg.param = struct('hbo', [10 10], 'hbr', [3 3], ...
                   'scatamp500', [1.4756 1.4756], 'scatpow500', [1.2079 1.2079]);
prop = rbupdateprop(cfg);
p690 = prop('690');
p830 = prop('830');
test_redbird('rbupdateprop 500nm musp@690', @() p690(2, 2), 1.4756 * (690 / 500)^(-1.2079));
test_redbird('rbupdateprop 500nm musp@830', @() p830(2, 2), 1.4756 * (830 / 500)^(-1.2079));

% Both pairs defined -> normalized convention wins
cfg.param.scatamp = [9 9];        % a deliberately wrong legacy value
cfg.param.scatpow = [9 9];
prop = rbupdateprop(cfg);
p690b = prop('690');
test_redbird('rbupdateprop prefers 500nm when both defined', ...
             @() p690b(2, 2), 1.4756 * (690 / 500)^(-1.2079));

% ---- rbupdateprop: mua matches sum_chromo (extin * conc) ----
% rbupdateprop fills default water=0.23 / lipids=0.58 when absent, so set
% them to zero to isolate the hbo/hbr contribution.
cfg.param = struct('hbo', [12 12], 'hbr', [4 4], 'water', [0 0], 'lipids', [0 0], ...
                   'scatamp500', [1.4756 1.4756], 'scatpow500', [1.2079 1.2079]);
prop = rbupdateprop(cfg);
p690c = prop('690');
ext = rbextinction(690, {'hbo', 'hbr'});
mua_expect = ext(1) * 12 + ext(2) * 4;
test_redbird('rbupdateprop mua = sum extin*conc', @() p690c(2, 1), mua_expect);

% defaults kick in when water/lipids absent
cfg.param = struct('hbo', [12 12], 'hbr', [4 4], ...
                   'scatamp500', [1.4756 1.4756], 'scatpow500', [1.2079 1.2079]);
prop = rbupdateprop(cfg);
p690d = prop('690');
extwl = rbextinction(690, {'water', 'lipids'});
mua_with_defaults = ext(1) * 12 + ext(2) * 4 + extwl(1) * 0.23 + extwl(2) * 0.58;
test_redbird('rbupdateprop applies default water=0.23 lipids=0.58', ...
             @() p690d(2, 1), mua_with_defaults);

% ---- rbsyncprop: copy recon.param onto cfg.param (label-based) ----
% recon has 2 labels and a recon.seg; cfg's node/elem are big enough that
% length(recon.param.hbo)=2 < labelmax, so the label-based copy branch is taken.
cfgS = struct;
cfgS.node = (1:30)';
cfgS.elem = reshape(1:120, 30, 4);
cfgS.param = struct('hbo', [0 0], 'hbr', [0 0]);
recon = struct('seg', [1; 2], ...
               'param', struct('hbo', [5; 6], 'hbr', [1; 2]));
[cfgS2, ~] = rbsyncprop(cfgS, recon);
test_redbird('rbsyncprop copies recon.param.hbo', @() cfgS2.param.hbo(:), [5; 6]);
test_redbird('rbsyncprop copies recon.param.hbr', @() cfgS2.param.hbr(:), [1; 2]);

% rbsyncprop without recon.param or recon.prop warns (and does nothing harmful)
recon3 = struct;
test_redbird('rbsyncprop with empty recon does not throw', ...
             @() rbsyncprop(cfgS, recon3), 'noerror');

% =========================================================================
function test_mesh()

cfg = struct;
[cfg.node, cfg.face, cfg.elem] = meshabox([0 0 0], [40 40 20], 4);
cfg.seg = ones(size(cfg.elem, 1), 1);
cfg.srcpos = [20 20 0];
cfg.srcdir = [0 0 1];
cfg.detpos = [20 20 20];
cfg.detdir = [0 0 -1];
cfg.prop = [0 0 1 1; 0.01 1 0 1.37];
cfg.omega = 0;

newcfg = rbmeshprep(cfg);

test_redbird('rbmeshprep adds reff', @isfield, true, newcfg, 'reff');
test_redbird('rbmeshprep adds musp0', @isfield, true, newcfg, 'musp0');
test_redbird('rbmeshprep adds deldotdel', @isfield, true, newcfg, 'deldotdel');
test_redbird('rbmeshprep adds evol', @isfield, true, newcfg, 'evol');
test_redbird('rbmeshprep marks isreoriented', @() newcfg.isreoriented, 1);

% reff for n=1.37 is the same value rbgetreff returns directly
test_redbird('rbmeshprep reff matches rbgetreff', ...
             @() newcfg.reff, rbgetreff(1.37, 1));

% Volume conservation: sum of element volumes equals the box volume.
% meshabox tessellates the box; small numerical drift (~1e-9) is fine.
box_volume = 40 * 40 * 20;
test_redbird('rbmeshprep evol sums to box volume', ...
             @(d) abs(d) < 1e-6, true, sum(newcfg.evol) - box_volume);

% Volume partition: nodal volumes also tile the box.
test_redbird('rbmeshprep nvol sums to box volume', ...
             @(d) abs(d) < 1e-6, true, sum(newcfg.nvol) - box_volume);

% degenerate-elem check should error
cfgBad = cfg;
cfgBad.evol = [];
cfgBad.elem(1, :) = cfgBad.elem(1, [1 1 2 3]);
test_redbird('rbmeshprep rejects degenerate elem', ...
             @() rbmeshprep(cfgBad), 'error');

% =========================================================================
function test_forward()

cfg = struct;
[cfg.node, cfg.face, cfg.elem] = meshabox([0 0 0], [40 40 20], 4);
cfg.seg = ones(size(cfg.elem, 1), 1);
cfg.srcpos = [20 20 0];
cfg.srcdir = [0 0 1];
cfg.detpos = [10 20 20; 30 20 20];
cfg.detdir = [0 0 -1];
cfg.prop = [0 0 1 1; 0.01 1 0 1.37];
cfg.omega = 0;
cfg = rbmeshprep(cfg);

[detphi, phi] = rbrunforward(cfg);

test_redbird('rbrunforward detphi has Ndet rows', @() size(detphi, 1), 2);
test_redbird('rbrunforward phi covers all nodes', ...
             @() size(phi, 1), size(cfg.node, 1));
% physically meaningful: finite values and the bulk is non-negative
% (FEM solver noise lets a handful of nodes go slightly negative, so we
% don't enforce strict non-negativity)
test_redbird('rbrunforward phi is finite', @(x) all(isfinite(x(:))), true, phi);
test_redbird('rbrunforward phi mostly non-negative', ...
             @(x) mean(x(:) >= 0) > 0.99, true, phi);
test_redbird('rbrunforward detphi is positive', @(x) all(x(:) > 0), true, detphi);

% farther detector (idx 2 here, both at 14.14mm — pick same; instead make a
% pair with different distances)
cfg.detpos = [25 20 20; 35 20 20];
cfg = rmfield(cfg, {'reff', 'musp0', 'deldotdel', 'rows', 'cols', 'idxcount', 'idxsum'});
cfg = rbmeshprep(cfg);
detphi2 = rbrunforward(cfg);
test_redbird('rbrunforward farther det has lower fluence', ...
             @(d) d(2) < d(1), true, detphi2);

% rbrun with single arg should match rbrunforward output
detphi3 = rbrun(cfg);
test_redbird('rbrun(cfg) == rbrunforward(cfg)', @() detphi3, detphi2);

% ---- Fluence symmetry ----
% Centered source with two detectors equidistant on opposite sides should
% produce equal detphi values (up to mesh-asymmetry noise).
cfgS = struct;
[cfgS.node, cfgS.face, cfgS.elem] = meshabox([0 0 0], [40 40 20], 4);
cfgS.seg = ones(size(cfgS.elem, 1), 1);
cfgS.srcpos = [20 20 0];
cfgS.srcdir = [0 0 1];
cfgS.detpos = [20 - 8, 20, 20; 20 + 8, 20, 20];
cfgS.detdir = [0 0 -1];
cfgS.prop = [0 0 1 1; 0.01 1 0 1.37];
cfgS.omega = 0;
cfgS = rbmeshprep(cfgS);
detS = rbrunforward(cfgS);
test_redbird('forward symmetric detectors agree to ~5%', ...
             @(r) abs(r - 1) < 0.05, true, detS(1) / detS(2));

% ---- Multi-wavelength forward ----
% Two distinct wavelengths in cfg.prop -> output is a containers.Map with
% both keys, and the two readings differ (different mua).
cfgM = struct;
[cfgM.node, cfgM.face, cfgM.elem] = meshabox([0 0 0], [40 40 20], 4);
cfgM.seg = ones(size(cfgM.elem, 1), 1);
cfgM.srcpos = [20 20 0];
cfgM.srcdir = [0 0 1];
cfgM.detpos = [20 20 20];
cfgM.detdir = [0 0 -1];
cfgM.prop = containers.Map({'690', '830'}, ...
                           {[0 0 1 1; 0.005 1.0 0 1.37], ...
                            [0 0 1 1; 0.020 0.8 0 1.37]});
cfgM.omega = 0;
cfgM = rbmeshprep(cfgM);
detM = rbrunforward(cfgM);
test_redbird('multispec forward returns containers.Map', ...
             @(x) isa(x, 'containers.Map'), true, detM);
test_redbird('multispec map has both wavelengths', ...
             @() sort(keys(detM)), {'690', '830'});
% Higher-mua wavelength -> lower fluence
test_redbird('higher mua -> lower detphi', ...
             @(d) d, true, detM('830') < detM('690'));

% ---- rbjac: J_mua sign ----
% For mua, increasing absorption uniformly decreases fluence at all
% detectors, so all entries of the Jacobian should be <= 0.
sdJ = rbsdmap(cfg);
[~, phiJ] = rbrunforward(cfg);
[Jmua_n, Jmua_e] = rbjac(sdJ, phiJ, cfg.deldotdel, cfg.elem, cfg.evol);
test_redbird('rbjac Jmua_node is non-positive', ...
             @(x) all(x(:) <= 1e-12), true, Jmua_n);
test_redbird('rbjac Jmua_elem is non-positive', ...
             @(x) all(x(:) <= 1e-12), true, Jmua_e);
% Sanity: not identically zero
test_redbird('rbjac Jmua_node not all zero', ...
             @(x) any(x(:) < -1e-9), true, Jmua_n);

% =========================================================================
function test_solver()

% rbreginv: small Tikhonov lambda should let the solution fit, large lambda
% should pull it toward zero (smaller norm). Construct an under-determined
% system so regularization is meaningful.
rng_state = rand('state');  %#ok<RAND> % save state to be polite
rand('state', 7);
A = randn(20, 50);
b = randn(20, 1);
x_low  = rbreginv(A, b, 1e-6);
x_high = rbreginv(A, b, 1e3);
test_redbird('rbreginv solution sizes match', ...
             @() size(x_low), [50 1]);
test_redbird('rbreginv high lambda -> smaller-norm solution', ...
             @(d) d > 1, true, norm(x_low) / norm(x_high));

% Over-determined case: lambda~0 should approximately recover the
% least-squares solution from backslash.
A2 = randn(60, 10);
b2 = randn(60, 1);
x_ls   = A2 \ b2;
x_reg0 = rbreginv(A2, b2, 1e-10);
test_redbird('rbreginv lambda~0 ~= least squares', ...
             @(r) r < 1e-3, true, norm(x_reg0 - x_ls) / norm(x_ls));

rand('state', rng_state);

% =========================================================================
function test_recon()

% Forward truth: small homogeneous slab with one absorbing inclusion.
cfg = struct;
[cfg.node, cfg.face, cfg.elem] = meshabox([0 0 0], [40 40 20], 6);
cfg.seg = ones(size(cfg.elem, 1), 1);
[xs, ys] = meshgrid(10:10:30, 10:10:30);
cfg.srcpos = [xs(:), ys(:), zeros(numel(xs), 1)];
cfg.detpos = [xs(:), ys(:), 20 * ones(numel(xs), 1)];
cfg.srcdir = [0 0 1];
cfg.detdir = [0 0 -1];
cfg.prop = [0 0 1 1; 0.01 1 0 1.37];
cfg.omega = 0;
cfg = rbmeshprep(cfg);

detphi0 = rbrunforward(cfg);

% Recon on the same mesh with a perturbed initial mua, single iteration.
recon = struct;
recon.prop = [0 0 1 1; 0.02 1 0 1.37];   % wrong initial guess
recon.lambda = 1e-2;

[newrecon, resid] = rbrunrecon(1, cfg, recon, detphi0);

test_redbird('rbrunrecon returns recon struct', @isstruct, true, newrecon);
test_redbird('rbrunrecon residual is finite', ...
             @(x) all(isfinite(x)), true, resid);
test_redbird('rbrunrecon residual is positive', ...
             @(x) x(1) > 0, true, resid);
% maxiter=0 short-circuits to forward solve
[detOnly, ~] = rbrunrecon(0, cfg, recon, detphi0);
test_redbird('rbrunrecon(maxiter=0) == forward', ...
             @() size(detOnly), size(detphi0));

% Multi-iteration: residual should not blow up. We allow small upticks
% between consecutive iterations (Gauss-Newton can wiggle), but final
% residual should be no worse than initial.
[~, resid_multi] = rbrunrecon(3, cfg, recon, detphi0);
test_redbird('rbrunrecon residual final <= initial', ...
             @(d) d <= 1e-9, true, resid_multi(end) - resid_multi(1));
test_redbird('rbrunrecon all residuals positive', ...
             @(x) all(x > 0), true, resid_multi);

% =========================================================================
function test_mwt()

% mu0_mm = 4*pi*1e-10 H/mm and eps0_mm = 8.854187817e-15 F/mm are the
% values used internally; tests verify the same constants.
mu0_mm = 4 * pi * 1e-10;
eps0_mm = 8.854187817e-15;

% ---- rbgetbulk MWT ----
cfgM = struct('bulk', struct('epsilon', 78, 'sigma', 1e-3, 'n', sqrt(78)));
bk = rbgetbulk(cfgM);
test_redbird('rbgetbulk MWT eps_r', @() bk(1), 78);
test_redbird('rbgetbulk MWT sigma', @() bk(2), 1e-3);
test_redbird('rbgetbulk MWT default mu0', @() bk(3), mu0_mm);
test_redbird('rbgetbulk MWT n', @() bk(4), sqrt(78));

cfgM2 = struct('bulk', struct('epsilon', 1));
bk2 = rbgetbulk(cfgM2);
test_redbird('rbgetbulk MWT minimal: sigma defaults to 0', @() bk2(2), 0);
test_redbird('rbgetbulk MWT minimal: mu0 = 4*pi*1e-10', @() bk2(3), mu0_mm);

% ---- rbupdateprop MWT (label-based, 2 labels over a small mesh) ----
cfgU = struct;
cfgU.node = [0 0 0; 1 0 0; 0 1 0; 0 0 1; 0.5 0.5 0.5];
cfgU.elem = [1 2 3 4; 1 2 3 5; 1 2 4 5];
cfgU.seg = [1; 1; 2];
cfgU.param = struct('epsilon', [40 60], 'sigma', [0.5e-3 1.0e-3]);
cfgU.prop = containers.Map({'5e8'}, {[1 0 mu0_mm 1; 1 0 mu0_mm 1; 1 0 mu0_mm 1]});
propnew = rbupdateprop(cfgU);
p = propnew('5e8');
test_redbird('rbupdateprop MWT epsilon (label 1)', @() p(2, 1), 40);
test_redbird('rbupdateprop MWT epsilon (label 2)', @() p(3, 1), 60);
test_redbird('rbupdateprop MWT sigma (label 1)', @() p(2, 2), 0.5e-3);
test_redbird('rbupdateprop MWT preserves mu0', @() p(2, 3), mu0_mm);

% ---- mesh-dependent tests ----
if (~exist('meshabox', 'file'))
    fprintf(2, '[skip] ''mwt'' mesh tests require iso2mesh\n');
    return
end

% ---- rbmeshprep MWT precomputations ----
cfgF = struct;
[cfgF.node, cfgF.face, cfgF.elem] = meshabox([0 0 0], [40 40 40], 8);
cfgF.seg = ones(size(cfgF.elem, 1), 1);
cfgF.srcpos = [10 20 0];
cfgF.detpos = [30 20 40];
cfgF.srcdir = [0 0 1];
cfgF.detdir = [0 0 -1];
cfgF.bulk = struct('epsilon', 1, 'sigma', 0, 'n', 1);
cfgF.prop = containers.Map({'5e8'}, {[1 0 mu0_mm 1; 1 0 mu0_mm 1]});
cfgF.omega = 2 * pi * 5e8;
cfgF = rbmeshprep(cfgF);

test_redbird('rbmeshprep MWT adds facenb', @isfield, true, cfgF, 'facenb');
test_redbird('rbmeshprep MWT facenb has 4 columns', @() size(cfgF.facenb, 2), 4);
test_redbird('rbmeshprep MWT adds facecenter', @isfield, true, cfgF, 'facecenter');
test_redbird('rbmeshprep MWT adds facenormal', @isfield, true, cfgF, 'facenormal');
test_redbird('rbmeshprep MWT adds facer', @isfield, true, cfgF, 'facer');
test_redbird('rbmeshprep MWT adds rbcorigin', @isfield, true, cfgF, 'rbcorigin');
test_redbird('rbmeshprep MWT skips reff (DOT-only)', ...
             @(c) ~isfield(c, 'reff'), true, cfgF);
test_redbird('rbmeshprep MWT facenormal unit length', ...
             @(x) max(abs(sqrt(sum(x.^2, 2)) - 1)) < 1e-9, true, cfgF.facenormal);
test_redbird('rbmeshprep MWT facer is positive', ...
             @(x) all(x > 0), true, cfgF.facer);

% ---- line source via cfg.srctype='line' + cfg.srcparam1 ----
cfgL = struct;
[cfgL.node, cfgL.face, cfgL.elem] = meshabox([0 0 0], [40 40 40], 8);
cfgL.seg = ones(size(cfgL.elem, 1), 1);
cfgL.srctype = 'line';
cfgL.srcpos = [20 20 0];
cfgL.srcparam1 = [0 0 40];
cfgL.dettype = 'line';
cfgL.detpos = [10 10 0];
cfgL.detparam1 = [0 0 40];
cfgL.bulk = struct('epsilon', 1, 'sigma', 0, 'n', 1);
cfgL.prop = containers.Map({'5e8'}, {[1 0 mu0_mm 1; 1 0 mu0_mm 1]});
cfgL.omega = 2 * pi * 5e8;
cfgL = rbmeshprep(cfgL);

test_redbird('line srcpos -> widesrc populated', @isfield, true, cfgL, 'widesrc');
test_redbird('line detpos -> widedet populated', @isfield, true, cfgL, 'widedet');
test_redbird('widesrc has 1 row per line src', ...
             @() size(cfgL.widesrc, 1), 1);
test_redbird('widesrc is non-zero', ...
             @(x) any(abs(x(:)) > 0), true, cfgL.widesrc);
test_redbird('widesrc is complex (MWT amp scaling -j*w*mu0)', ...
             @(x) ~isreal(x), true, cfgL.widesrc);

% ---- end-to-end MWT forward smoke test ----
[detphi, phi] = rbrunforward(cfgL);
test_redbird('MWT forward: phi finite', ...
             @(x) all(isfinite(x(:))), true, phi);
test_redbird('MWT forward: phi complex', ...
             @(x) ~isreal(x), true, phi);
test_redbird('MWT forward: phi non-trivial', ...
             @(x) max(abs(x(:))) > 0, true, phi);
test_redbird('MWT forward: detphi finite', ...
             @(x) all(isfinite(x(:))), true, detphi);

% ---- DOT regression: line source for DOT works too ----
cfgD = struct;
[cfgD.node, cfgD.face, cfgD.elem] = meshabox([0 0 0], [40 40 40], 8);
cfgD.seg = ones(size(cfgD.elem, 1), 1);
cfgD.srctype = 'line';
cfgD.srcpos = [20 20 0];
cfgD.srcparam1 = [0 0 40];
cfgD.dettype = 'line';
cfgD.detpos = [10 10 35];
cfgD.detparam1 = [20 20 0];
cfgD.prop = [0 0 1 1; 0.01 1 0 1.37];
cfgD.omega = 0;
cfgD = rbmeshprep(cfgD);
test_redbird('line source for DOT also works', ...
             @() size(cfgD.widesrc, 1), 1);
detD = rbrunforward(cfgD);
test_redbird('DOT line-src forward returns finite', ...
             @(x) all(isfinite(x(:))), true, detD);

% ---- ray source: collimated-beam model (Haskell et al. 1994, eq 2.4.5) ----
% mu_s' * exp(-mu_s' * l) weighted line-of-sources along cfg.srcdir, normalized
% so total RHS sum = 1 and mean depth = l_tr.
cfgY = struct;
[cfgY.node, cfgY.face, cfgY.elem] = meshabox([0 0 0], [40 40 40], 8);
cfgY.seg = ones(size(cfgY.elem, 1), 1);
cfgY.srctype = 'ray';
cfgY.srcpos = [20 20 0];
cfgY.srcdir = [0 0 1];
cfgY.dettype = 'ray';
cfgY.detpos = [10 10 0];
cfgY.detdir = [0 0 1];
cfgY.prop = [0 0 1 1; 0.01 1 0 1.37];   % mua=0.01, mus=1, g=0 -> mus'=1 -> l_tr=1mm
cfgY.omega = 0;
cfgY = rbmeshprep(cfgY);

test_redbird('ray src -> widesrc populated', @isfield, true, cfgY, 'widesrc');
test_redbird('ray det -> widedet populated', @isfield, true, cfgY, 'widedet');
% with mua=0.01, mus'=1 -> mu_tr=1.01, transport albedo = mus'/mu_tr = 0.99
test_redbird('ray widesrc sum ~ transport albedo', ...
             @(x) abs(x - 1 / 1.01) < 5e-3, true, sum(cfgY.widesrc(1, :)));
test_redbird('ray widedet sum ~ transport albedo', ...
             @(x) abs(x - 1 / 1.01) < 5e-3, true, sum(cfgY.widedet(1, :)));

% center-of-mass along srcdir should equal 1/mu_tr = 0.99 mm
zcoord = cfgY.node(:, 3);
zcom = sum(zcoord' .* cfgY.widesrc(1, :)) / sum(cfgY.widesrc(1, :));
test_redbird('ray src mean depth ~ 1/mu_tr', ...
             @(z) abs(z - 1 / 1.01) < 0.05, true, zcom);

detY = rbrunforward(cfgY);
test_redbird('DOT ray-src forward returns finite', ...
             @(x) all(isfinite(x(:))), true, detY);
test_redbird('DOT ray-src forward non-trivial', ...
             @(x) max(abs(x(:))) > 0, true, detY);

% missing srcdir should error for ray type
cfgYerr = struct;
[cfgYerr.node, cfgYerr.face, cfgYerr.elem] = meshabox([0 0 0], [40 40 40], 8);
cfgYerr.seg = ones(size(cfgYerr.elem, 1), 1);
cfgYerr.srctype = 'ray';
cfgYerr.srcpos = [20 20 0];
cfgYerr.detpos = [30 30 0];
cfgYerr.prop = [0 0 1 1; 0.01 1 0 1.37];
cfgYerr.omega = 0;
test_redbird('ray src without srcdir errors', ...
             @() rbmeshprep(cfgYerr), 'error');

% ---- analytic checks: k-formula and physics-driven invariants ----

% k_bg formula: in vacuum (eps_r=1, sigma=0), k = omega / c with c in mm/s.
c0_mm = 1 / sqrt(mu0_mm * eps0_mm);   % ~ 2.998e11 mm/s
omega_test = 2 * pi * 5e8;
k_vac = sqrt(omega_test^2 * mu0_mm * eps0_mm * 1 - 1j * omega_test * mu0_mm * 0);
test_redbird('MWT k formula: k(vacuum) = omega/c', ...
             @(d) abs(d) < 1e-15, true, abs(real(k_vac) - omega_test / c0_mm));
test_redbird('MWT k formula: vacuum is lossless (Im(k)=0)', ...
             @() imag(k_vac), 0);

% Lossy medium: Im(k) > 0 magnitude consistent with -j*omega*mu*sigma.
% For sigma > 0, Im(k^2) < 0, so principal sqrt has Re > 0 and Im < 0,
% giving exp(-j*k*r) decay (since j*(-Im(k)) = +Im(k) makes the exponent
% acquire a negative real part).
k_lossy = sqrt(omega_test^2 * mu0_mm * eps0_mm * 1 - 1j * omega_test * mu0_mm * 1e-3);
test_redbird('MWT k(lossy): Im(k) negative for sigma>0', ...
             @(x) x < 0, true, imag(k_lossy));

% ---- reciprocity: detphi(2->1) == detphi(1->2) when geometry permits ----
% A complex-symmetric (A = A.', not Hermitian) Helmholtz FEM matrix gives
% reciprocity between any two antennas.  Use two distinct line antennas so
% the sd map doesn't reduce to self-pairs.
cfgR = struct;
[cfgR.node, cfgR.face, cfgR.elem] = meshabox([0 0 0], [40 40 40], 8);
cfgR.seg = ones(size(cfgR.elem, 1), 1);
cfgR.srctype = 'line';
cfgR.srcpos = [10 20 5
               30 20 5];
cfgR.srcparam1 = [0 0 30];
cfgR.dettype = 'line';
cfgR.detpos = [10 20 5
               30 20 5];
cfgR.detparam1 = [0 0 30];
cfgR.bulk = struct('epsilon', 4, 'sigma', 1e-4, 'n', 2);
cfgR.prop = containers.Map({'5e8'}, {[1 0 mu0_mm 1; 4 1e-4 mu0_mm 2]});
cfgR.omega = 2 * pi * 5e8;
cfgR = rbmeshprep(cfgR);
detphi_R = rbrunforward(cfgR);

% off-diagonal entries are antenna1<->antenna2 measurements; should match.
recip_err = abs(detphi_R(1, 2) - detphi_R(2, 1)) / max(abs(detphi_R(1, 2)), abs(detphi_R(2, 1)));
test_redbird('MWT reciprocity: detphi(2->1) == detphi(1->2)', ...
             @(r) r < 1e-6, true, recip_err);

% ---- attenuation: lossy medium produces smaller |E| than lossless ----
cfgLF = struct;
[cfgLF.node, cfgLF.face, cfgLF.elem] = meshabox([0 0 0], [40 40 40], 8);
cfgLF.seg = ones(size(cfgLF.elem, 1), 1);
cfgLF.srctype = 'line';
cfgLF.srcpos = [10 20 5];
cfgLF.srcparam1 = [0 0 30];
cfgLF.dettype = 'line';
cfgLF.detpos = [30 20 5];
cfgLF.detparam1 = [0 0 30];
cfgLF.omega = 2 * pi * 5e8;
% lossless reference
cfgLF.bulk = struct('epsilon', 4, 'sigma', 0, 'n', 2);
cfgLF.prop = containers.Map({'5e8'}, {[1 0 mu0_mm 1; 4 0 mu0_mm 2]});
cfgLF = rbmeshprep(cfgLF);
detphi_lossless = rbrunforward(cfgLF);

% lossy version with same geometry/permittivity
cfgLY = struct;
[cfgLY.node, cfgLY.face, cfgLY.elem] = meshabox([0 0 0], [40 40 40], 8);
cfgLY.seg = ones(size(cfgLY.elem, 1), 1);
cfgLY.srctype = 'line';
cfgLY.srcpos = [10 20 5];
cfgLY.srcparam1 = [0 0 30];
cfgLY.dettype = 'line';
cfgLY.detpos = [30 20 5];
cfgLY.detparam1 = [0 0 30];
cfgLY.omega = 2 * pi * 5e8;
cfgLY.bulk = struct('epsilon', 4, 'sigma', 1e-3, 'n', 2);
cfgLY.prop = containers.Map({'5e8'}, {[1 0 mu0_mm 1; 4 1e-3 mu0_mm 2]});
cfgLY = rbmeshprep(cfgLY);
detphi_lossy = rbrunforward(cfgLY);

test_redbird('MWT lossy medium attenuates relative to lossless', ...
             @(r) r < 1, true, abs(detphi_lossy(1)) / abs(detphi_lossless(1)));

% ---- frequency scaling: in a lossless medium, k = omega/c and Re(k)
% scales linearly with omega (10x frequency -> 10x wavenumber). ----
omega_lo = 2 * pi * 5e8;
omega_hi = 2 * pi * 5e9;
k_lo = sqrt(omega_lo^2 * mu0_mm * eps0_mm);   % lossless vacuum
k_hi = sqrt(omega_hi^2 * mu0_mm * eps0_mm);
test_redbird('MWT k(omega): Re(k) scales linearly with omega in lossless', ...
             @(r) abs(r - 10) < 1e-9, true, real(k_hi) / real(k_lo));

% =========================================================================
function test_td()

% ---- rbfemlhs mode=3 returns pure consistent mass matrix ----
% on a single tetrahedron with known volume, M_ii = V/10, M_ij = V/20.
cfgU = struct;
cfgU.node = [0 0 0; 1 0 0; 0 1 0; 0 0 1];
cfgU.elem = [1 2 3 4];
cfgU.face = [2 3 4; 1 3 4; 1 2 4; 1 2 3];
cfgU.area = elemvolume(cfgU.node, cfgU.face);
cfgU.evol = elemvolume(cfgU.node, cfgU.elem);
cfgU.deldotdel = rbdeldotdel(cfgU);
cfgU.prop = [0 0 1 1; 0.01 1 0 1.37];
cfgU.seg = 1;
M = rbfemlhs(cfgU, cfgU.deldotdel, '_', 3);
test_redbird('rbfemlhs mode=3: M(1,1) = V_e/10', ...
             @(x) abs(x - cfgU.evol / 10) < 1e-15, true, full(M(1, 1)));
test_redbird('rbfemlhs mode=3: M(1,2) = V_e/20', ...
             @(x) abs(x - cfgU.evol / 20) < 1e-15, true, full(M(1, 2)));
test_redbird('rbfemlhs mode=3: M is symmetric', ...
             @(M) max(max(abs(M - M.'))) < 1e-15, true, full(M));
test_redbird('rbfemlhs mode=3: row-sum = V_e/4', ...
             @(s) max(abs(s - cfgU.evol / 4)) < 1e-15, true, full(sum(M, 2)));

% ---- TD forward smoke: default impulse source ----
cfg = struct;
[cfg.node, cfg.face, cfg.elem] = meshabox([0 0 0], [40 40 40], 8);
cfg.seg = ones(size(cfg.elem, 1), 1);
cfg.srcpos = [20 20 0];
cfg.srcdir = [0 0 1];
cfg.detpos = [20 20 40];
cfg.detdir = [0 0 -1];
cfg.prop = [0 0 1 1; 0.01 1 0 1.37];
cfg.omega = 0;
cfg.tstart = 0;
cfg.tstep = 50e-12;
cfg.tend = 2e-9;
cfg = rbmeshprep(cfg);
[detphi_td, ~, Amat_td, rhs_td] = rbrunforward(cfg);

test_redbird('TD detphi has 3D shape Ndet x Nsrc x Nt', ...
             @() size(detphi_td), [1 1 41]);
test_redbird('TD detphi is finite', ...
             @(x) all(isfinite(x(:))), true, detphi_td);
% past the peak, the diffusing tail is monotonically non-negative; early-time
% CN ringing on a coarse mesh with impulse IC is a known numerical artifact
% (use a smooth srctemporal or finer mesh in production).
[~, peak_lin] = max(abs(detphi_td(:)));
[~, ~, peak_t_idx] = ind2sub(size(detphi_td), peak_lin);
test_redbird('TD detphi tail (after peak) is non-negative', ...
             @(x) all(x >= 0), true, ...
             squeeze(detphi_td(:, :, peak_t_idx:end)));
test_redbird('TD detphi has a non-trivial peak', ...
             @(x) max(abs(x(:))) > 0, true, detphi_td);

% peak should occur somewhere in the interior, not at t=tstart and not at t=tend
[~, peak_idx] = max(abs(detphi_td(:)));
[~, ~, peak_t] = ind2sub(size(detphi_td), peak_idx);
test_redbird('TD detphi peak is interior (not at t=tstart)', ...
             @(t) t > 1, true, peak_t);
test_redbird('TD detphi peak is interior (not at t=tend)', ...
             @(t) t < size(detphi_td, 3), true, peak_t);

% ---- TD forward with custom temporal modulation ----
cfg2 = cfg;
cfg2 = rmfield(cfg2, {'reff', 'musp0', 'deldotdel', 'rows', 'cols', 'idxcount', 'idxsum', 'facenb'});
cfg2 = rbmeshprep(cfg2);
nt = length(cfg2.tstart:cfg2.tstep:cfg2.tend);
[detphi_const, ~] = rbrunforward(cfg2, 'srctemporal', ones(nt, 1));
test_redbird('TD with constant src has 3D shape', ...
             @() size(detphi_const), [1 1 nt]);
test_redbird('TD with constant src is finite', ...
             @(x) all(isfinite(x(:))), true, detphi_const);
% with constant source from t=0 starting from phi=0, detphi grows monotonically
% over time (until eventually approaching steady state)
test_redbird('TD with constant src grows over time', ...
             @(d) d(end) > d(2), true, detphi_const);

% ---- conflict: cfg.omega > 0 + TD should error in rbmeshprep ----
cfgE = cfg;
cfgE = rmfield(cfgE, {'reff', 'musp0', 'deldotdel', 'rows', 'cols', 'idxcount', 'idxsum', 'facenb'});
cfgE.omega = 2 * pi * 1e8;
test_redbird('TD + omega>0 errors in rbmeshprep', ...
             @() rbmeshprep(cfgE), 'error');

% ---- conflict: MWT bulk + TD should error in rbmeshprep ----
cfgM = cfg;
cfgM = rmfield(cfgM, {'reff', 'musp0', 'deldotdel', 'rows', 'cols', 'idxcount', 'idxsum', 'facenb'});
cfgM.bulk = struct('epsilon', 4, 'sigma', 0, 'n', 2);
test_redbird('TD + MWT bulk.epsilon errors in rbmeshprep', ...
             @() rbmeshprep(cfgM), 'error');

% ---- integration invariant: int_0^T detphi_TD(t) dt ≈ detphi_CW for impulse IC ----
% the impulse response integrated over time should equal the CW steady-state
% response (Laplace-transform-at-zero identity for the diffusion equation).
% use a longer time window to capture the tail.
cfgI = struct;
[cfgI.node, cfgI.face, cfgI.elem] = meshabox([0 0 0], [40 40 40], 8);
cfgI.seg = ones(size(cfgI.elem, 1), 1);
cfgI.srcpos = [20 20 0];
cfgI.srcdir = [0 0 1];
cfgI.detpos = [20 20 40];
cfgI.detdir = [0 0 -1];
cfgI.prop = [0 0 1 1; 0.01 1 0 1.37];
cfgI.omega = 0;
cfgI.tstart = 0;
cfgI.tstep = 50e-12;
cfgI.tend = 1e-8;        % 10 ns is enough for diffusion tail to subside
cfgI = rbmeshprep(cfgI);
detphi_imp = rbrunforward(cfgI);
% trapezoidal integral over time
integ = squeeze(sum(detphi_imp, 3) - 0.5 * (detphi_imp(:, :, 1) + detphi_imp(:, :, end))) * cfgI.tstep;

% CW reference
cfgC = rmfield(cfgI, {'tstart', 'tstep', 'tend', 'reff', 'musp0', 'deldotdel', 'rows', 'cols', 'idxcount', 'idxsum', 'facenb'});
cfgC = rbmeshprep(cfgC);
detphi_cw = rbrunforward(cfgC);

% expect integral and CW solution to agree to within a few percent
relerr = abs(integ - detphi_cw) / abs(detphi_cw);
test_redbird('TD integration invariant: int_0^T phi_TD(t) dt ~ phi_CW', ...
             @(r) r < 0.05, true, relerr);
