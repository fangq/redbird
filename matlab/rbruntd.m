function [detphi, phi, Amat, rhs] = rbruntd(cfg, varargin)
%
% [detphi, phi, Amat, rhs]=rbruntd(cfg)
%   or
% [detphi, phi, Amat, rhs]=rbruntd(cfg,'param1',value1,...)
%
% Time-domain DOT forward solver using an implicit Crank-Nicolson scheme
% for the time-dependent diffusion equation:
%
%   -div(D grad Phi) + mua*Phi + (1/c) dPhi/dt = S(r,t)
%
% engaged automatically by rbrunforward when cfg.tstart, cfg.tstep and
% cfg.tend are all defined.
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     cfg: the redbird simulation data structure with the standard DOT
%          fields plus the time-domain triplet:
%          cfg.tstart - start time (s)
%          cfg.tstep  - time step Delta_t (s)
%          cfg.tend   - end time (s)
%     param/value pairs:
%          'theta':       theta-method weight (default 0.5 = Crank-Nicolson;
%                         set to 1.0 for backward Euler in stiff regimes)
%          'srctemporal': temporal modulation of the source. defaults to an
%                         impulse at t=tstart (TPSF interpretation; the
%                         initial condition is set to c*M\\S_spatial). Can
%                         be a length-Nt vector, a function handle f(t),
%                         or omitted for the impulse default.
%          'phi0':        initial condition Nn x Nsrc; defaults to zeros
%                         (or the impulse IC when srctemporal is empty).
%          'tdsavevol':   1 to also return the full volumetric phi (Nn x
%                         Nsrc x Nt); 0 (default) returns only detphi.
%          'sd':          source-detector mapping (otherwise rbsdmap is
%                         called).
%
% output:
%     detphi: detector readings, Ndet x Nsrc x Nt array (single
%             wavelength) or a containers.Map keyed by wavelength.
%     phi:    volumetric forward solution, Nn x Nsrc x Nt (only when
%             tdsavevol is true; otherwise [] or empty Map entries).
%     Amat:   the time-step LHS A_TD = M/(c*dt) + theta*A_cw, single matrix
%             or containers.Map keyed by wavelength.
%     rhs:    the spatial source/detector right-hand-side from rbfemrhs.
%
% license:
%     GPL version 3, see LICENSE_GPLv3.txt files for details
%
% -- this function is part of Redbird-m toolbox
%

opt = varargin2struct(varargin{:});
theta = jsonopt('theta', 0.5, opt);
tdsavevol = jsonopt('tdsavevol', 0, opt);
phi0 = jsonopt('phi0', [], opt);
srctemporal = jsonopt('srctemporal', [], opt);

if (~isfield(cfg, 'tstart') || ~isfield(cfg, 'tstep') || ~isfield(cfg, 'tend'))
    error('rbruntd requires cfg.tstart, cfg.tstep, and cfg.tend');
end
if (cfg.tend <= cfg.tstart || cfg.tstep <= 0)
    error('cfg.tend must exceed cfg.tstart and cfg.tstep must be positive');
end

% time grid (inclusive endpoint)
tvec = cfg.tstart:cfg.tstep:cfg.tend;
nt = length(tvec);

% bulk refractive index for the 1/c prefactor on the time derivative
bk = rbgetbulk(cfg);
if (isa(bk, 'containers.Map'))
    kk = bk.keys;
    bk = bk(kk{1});
end
nref_bulk = bk(4);
c_mm_per_s = 299792458000 / nref_bulk;

% wavelength keys (use the same '_' placeholder pattern as rbrunforward)
wavelengths = {'_'};
if (isa(cfg.prop, 'containers.Map'))
    wavelengths = cfg.prop.keys;
end

% sd map
if isfield(opt, 'sd') && ~isempty(opt.sd)
    sd = opt.sd;
else
    sd = rbsdmap(cfg);
end
if (~isa(sd, 'containers.Map'))
    sd = containers.Map(wavelengths, {sd});
end

if (~isfield(cfg, 'deldotdel'))
    cfg.deldotdel = rbdeldotdel(cfg);
end

nn = size(cfg.node, 1);

% identify source/detector column ranges of rbfemrhs output. the rhs is
% Nn x (Nsrc + Ndet) with sources first, detectors last (rbfemrhs line 85).
srcnum = 0;
if (isfield(cfg, 'srcpos') && ~isempty(cfg.srcpos))
    srcnum = srcnum + size(cfg.srcpos, 1);
end
if (isfield(cfg, 'widesrc') && ~isempty(cfg.widesrc))
    srcnum = srcnum + size(cfg.widesrc, 1);
end
detnum = 0;
if (isfield(cfg, 'detpos') && ~isempty(cfg.detpos))
    detnum = detnum + size(cfg.detpos, 1);
end
if (isfield(cfg, 'widedet') && ~isempty(cfg.widedet))
    detnum = detnum + size(cfg.widedet, 1);
end

Amat = containers.Map();
phi = containers.Map();
detphi = containers.Map();
rhs = [];

% the consistent mass matrix M_ij = <phi_i, phi_j> is wavelength-independent
% (geometry only: cfg.elem and cfg.evol), so build it once outside the loop.
M = rbfemlhs(cfg, cfg.deldotdel, wavelengths{1}, 3);

for waveid = wavelengths
    wv = waveid{1};

    % spatial CW operator A_cw = D*K + mua*M + Robin BC (per wavelength)
    A_cw = rbfemlhs(cfg, cfg.deldotdel, wv, 2);

    % Crank-Nicolson time-step operators (constant Delta_t -> reuse factorization)
    M_T  = M / (c_mm_per_s * cfg.tstep);
    A_TD = M_T + theta * A_cw;
    B_TD = M_T - (1 - theta) * A_cw;

    Amat(wv) = A_TD;

    % factor LHS once per wavelength (constant Delta_t -> constant LHS).
    % use decomposition() in MATLAB for the cleanest API; fall back to a
    % pre-computed sparse LU in Octave.
    use_decomp = exist('decomposition', 'builtin') || exist('decomposition', 'file');
    if (use_decomp)
        dA_TD = decomposition(A_TD);
    else
        [LL, UU, PP, QQ] = lu(A_TD);
    end

    % spatial source/detector vectors (Nn x (srcnum+detnum))
    rhs_spatial = rbfemrhs(cfg, sd, wv, 1);
    if (size(rhs_spatial, 2) ~= srcnum + detnum)
        srcnum = size(rhs_spatial, 2) - detnum;
    end
    S_spatial = full(rhs_spatial(:, 1:srcnum));
    D_spatial = full(rhs_spatial(:, srcnum + 1:srcnum + detnum));

    rhs = rhs_spatial;

    % temporal modulation
    if (isempty(srctemporal))
        srct = zeros(nt, 1);   % impulse: handled via IC, no driving source
        is_impulse = true;
    elseif (isa(srctemporal, 'function_handle'))
        srct = srctemporal(tvec(:));
        is_impulse = false;
    else
        srct = srctemporal(:);
        if (length(srct) ~= nt)
            error('cfg.srctemporal must have %d entries (one per time step)', nt);
        end
        is_impulse = false;
    end

    % initial condition. for impulse default: Phi(t_start) = c * M\S_spatial
    % (the unique solution to the impulse-jump equation M*dPhi = c*S_spatial*dt
    % integrated across the delta).
    if (~isempty(phi0))
        phi_prev = phi0;
    elseif (is_impulse)
        phi_prev = c_mm_per_s * (M \ S_spatial);
    else
        phi_prev = zeros(nn, srcnum);
    end

    % output allocation
    if (tdsavevol)
        phi_t = zeros(nn, srcnum, nt);
        phi_t(:, :, 1) = phi_prev;
    end
    detphi_t = zeros(detnum, srcnum, nt);
    if (detnum > 0 && srcnum > 0)
        detphi_t(:, :, 1) = D_spatial' * phi_prev;
    end

    % time-stepping loop (n indexes the new time level)
    for nstep = 2:nt
        if (is_impulse)
            rhs_step = B_TD * phi_prev;
        else
            rhs_step = B_TD * phi_prev + 0.5 * (srct(nstep - 1) + srct(nstep)) * S_spatial;
        end
        if (use_decomp)
            phi_new = dA_TD \ rhs_step;
        else
            phi_new = QQ * (UU \ (LL \ (PP * rhs_step)));
        end
        if (tdsavevol)
            phi_t(:, :, nstep) = phi_new;
        end
        if (detnum > 0 && srcnum > 0)
            detphi_t(:, :, nstep) = D_spatial' * phi_new;
        end
        phi_prev = phi_new;
    end

    detphi(wv) = detphi_t;
    if (tdsavevol)
        phi(wv) = phi_t;
    end
end

% unwrap single-wavelength Map for ergonomic returns
if (length(wavelengths) == 1)
    Amat = Amat(wavelengths{1});
    detphi = detphi(wavelengths{1});
    if (tdsavevol)
        phi = phi(wavelengths{1});
    else
        phi = [];
    end
end
