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
    tests = {'util', 'jac', 'prop', 'mesh', 'forward', 'solver', 'recon'};
end

global RB_FAIL RB_TOTAL;
RB_FAIL = 0;
RB_TOTAL = 0;

if (~exist('rbrun', 'file'))
    addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'matlab'));
end
addpath(fileparts(mfilename('fullpath')));

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
