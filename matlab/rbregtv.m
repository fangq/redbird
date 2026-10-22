function [delta_mu, info] = rbregtv(J, r, alpha, varargin)
%
% [delta_mu, info] = rbregtv(J, r, alpha, 'param', value, ...)
%
% Total-variation regularized least-squares reconstruction via FISTA.
% Solves the composite-objective problem
%
%     min_x   0.5 * || J * x - r ||_2^2   +   alpha * TV_iso(x)
%
% where TV_iso(x) is the isotropic 3D total-variation seminorm
%
%     TV_iso(x) = sum_voxel sqrt((D_x x)^2 + (D_y x)^2 + (D_z x)^2)
%
% via Beck & Teboulle (2009)-style FISTA (Fast Iterative Shrinkage-
% Thresholding Algorithm).  Each outer FISTA step
%   - takes one gradient step on the smooth data-fit term, then
%   - applies the TV proximal operator (a 3D TV-denoising sub-problem
%     solved by ~10 inner iterations of Chambolle's (2004) projected-
%     gradient algorithm).
%
% Unlike rbreglsqr (minimum-L2-norm), TV-regularization biases the
% reconstruction toward piecewise-constant solutions, which much better
% preserves compact inclusions in DOT mua maps.  The cost is more
% iterations + an inner loop, and one tuning parameter (alpha) that
% trades reconstruction smoothness against measurement fit.
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     J:      4D array (Nx, Ny, Nz, Nsd) or function handle afun(x, flag);
%             same convention as rbreglsqr / rbjacop.
%     r:      measurement residual (Nsd-vector).
%     alpha:  TV regularization weight (>= 0).  Typical range:
%             1e-5 ... 1e-2 relative to ||J^T r||_inf.  Larger alpha ->
%             flatter / blockier recon.  Tune via an L-curve or by
%             matching the discrepancy principle.
%
%   param/value pairs (all optional):
%     'maxit'        : outer FISTA iters (default 100).
%     'innerit'      : Chambolle TV-prox inner iters per FISTA step
%                      (default 10).
%     'tol'          : relative-change stopping tolerance on ||x_{k+1} -
%                      x_k|| / ||x_k|| (default 1e-5).
%     'L'            : Lipschitz constant of J^T J (default: power-iter
%                      estimate; pass a known value to skip estimation).
%     'reshape3d'    : reshape result to (Nx, Ny, Nz) when J is 4D
%                      (default true).
%     'adjointtest'  : run one-shot <Jv, u> == <v, J' u> dot-product
%                      check (default true).
%     'verbose'      : print per-iteration progress (default 0).
%     'positive'     : project x >= 0 inside each FISTA prox step
%                      (default false).  Use this for mua perturbation
%                      recon where the physically meaningful delta_mua
%                      is non-negative (e.g. an absorbing inclusion in
%                      a uniform background).  Combined with TV, the
%                      positivity prior dramatically reduces sign-
%                      ambiguous noise-fitting at the array edges.
%
% output:
%     delta_mu: reconstructed voxel update (Nx, Ny, Nz) or flat Nv-vector.
%     info:     struct with fields
%                  .iter         outer FISTA iterations taken
%                  .L            Lipschitz constant used
%                  .obj          objective-value history
%                  .relchg       per-iter ||x_{k+1} - x_k|| / ||x_k||
%                  .adjoint_err  relative error in the dot-product test
%
% see also: rbreglsqr (L2 min-norm), rbjacop, rbreginv
%
% references:
%   * Beck & Teboulle, "A Fast Iterative Shrinkage-Thresholding
%     Algorithm for Linear Inverse Problems," SIAM J. Imaging Sci.,
%     2(1):183-202 (2009).
%   * Chambolle, "An Algorithm for Total Variation Minimization and
%     Applications," J. Math. Imaging Vis., 20(1-2):89-97 (2004).
%
% license:
%     GPL version 3, see LICENSE_GPLv3.txt files for details
%
% -- this function is part of Redbird toolbox
%

opt = varargin2struct(varargin{:});
maxit      = jsonopt('maxit',        100,   opt);
innerit    = jsonopt('innerit',      10,    opt);
tol        = jsonopt('tol',          1e-5,  opt);
L_in       = jsonopt('L',            [],    opt);
do_reshape = jsonopt('reshape3d',    true,  opt);
do_adjtst  = jsonopt('adjointtest',  true,  opt);
verbose    = jsonopt('verbose',      0,     opt);
nonneg     = jsonopt('positive',     false, opt);   % project x >= 0 each step

% --- wrap J ---
input_was_4d = (isnumeric(J) && ndims(J) == 4);
if (input_was_4d)
    sz4 = size(J);
    dims3 = sz4(1:3);
    Nv = prod(dims3);
    Nsd = sz4(4);
    afun = rbjacop(J);
elseif (isa(J, 'function_handle'))
    afun = J;
    if (isempty(r))
        error('rbregtv: residual r must be non-empty to size the operator');
    end
    Nsd = numel(r);
    probe = afun(zeros(Nsd, 1), 'transp');
    Nv = numel(probe);
    dims3 = jsonopt('dims3', [], opt);
    if (isempty(dims3))
        error(['rbregtv: when J is a function handle, you must pass ' ...
               '''dims3'', [Nx Ny Nz] so the TV operators know the volume shape.']);
    end
elseif (isstruct(J) && isfield(J, 'mua'))
    Jarr = J.mua;
    sz4 = size(Jarr);
    dims3 = sz4(1:3);
    Nv = prod(dims3);
    Nsd = sz4(4);
    afun = rbjacop(Jarr);
    input_was_4d = true;
else
    error('rbregtv: J must be a 4D array, a function handle, or a struct with field .mua');
end

% --- adjoint dot-product test ---
adjoint_err = NaN;
if (do_adjtst)
    v_probe = randn(Nv, 1);
    u_probe = randn(Nsd, 1);
    Jv  = afun(v_probe, 'notransp');
    Jtu = afun(u_probe, 'transp');
    lhs = u_probe(:).' * Jv(:);
    rhs = Jtu(:).' * v_probe(:);
    adjoint_err = abs(lhs - rhs) / max(abs(lhs), eps);
    if (adjoint_err > 1e-4)
        warning('rbregtv:adjointMismatch', ...
                'Adjoint dot-product test failed: rel err = %.3e', adjoint_err);
    end
end

% --- Lipschitz constant via power iteration on J^T J ---
if (isempty(L_in))
    L = power_iter_jtj(afun, Nv, 10);
    L = L * 1.05;     % small safety factor
else
    L = L_in;
end

% --- FISTA outer loop ---
r_col = r(:);
x_prev = zeros(Nv, 1);
y = zeros(Nv, 1);
t = 1.0;

obj_hist = zeros(maxit, 1);
relchg_hist = zeros(maxit, 1);

for k = 1:maxit
    % gradient of the data-fit term at y
    Jy = afun(y, 'notransp');
    grad = afun(Jy - r_col, 'transp');

    % gradient step
    z = y - grad / L;

    % proximal step: 3D TV denoising on z; optionally project x >= 0
    z3 = reshape(z, dims3);
    x_new3 = tv_prox_iso3d(z3, alpha / L, innerit);
    if (nonneg)
        x_new3 = max(x_new3, 0);
    end
    x_new = x_new3(:);

    % FISTA acceleration
    t_new = 0.5 * (1 + sqrt(1 + 4 * t * t));
    y = x_new + ((t - 1) / t_new) * (x_new - x_prev);

    % bookkeeping
    res = Jy - r_col;
    obj_hist(k) = 0.5 * (res(:).' * res(:)) + alpha * tv_iso3d_value(x_new3);
    relchg_hist(k) = norm(x_new - x_prev) / max(norm(x_prev), eps);

    if (verbose && (k <= 5 || mod(k, 10) == 0 || k == maxit))
        fprintf('  FISTA iter %3d:  obj=%.4e   relchg=%.3e\n', k, obj_hist(k), relchg_hist(k));
    end

    if (k > 1 && relchg_hist(k) < tol)
        break
    end

    x_prev = x_new;
    t = t_new;
end

niter = k;

info = struct('iter', niter, 'L', L, ...
              'obj', obj_hist(1:niter), 'relchg', relchg_hist(1:niter), ...
              'adjoint_err', adjoint_err);

if (input_was_4d && do_reshape)
    delta_mu = reshape(x_new, dims3);
else
    delta_mu = x_new;
end

end

% ---------------------------------------------------------------------
% helpers below: power-iteration Lipschitz estimate + 3D isotropic
% TV prox via Chambolle (2004) projected gradient + forward/back diff.
% ---------------------------------------------------------------------

function L = power_iter_jtj(afun, Nv, niter)
% Estimate ||J^T J||_2 via power iteration.  niter ~ 10 is enough for
% a usable upper bound; FISTA is robust to slight over-estimation of L
% (just slows convergence linearly).
v = randn(Nv, 1);
v = v / norm(v);
for k = 1:niter
    w = afun(v, 'notransp');
    v = afun(w, 'transp');
    nrm = norm(v);
    if (nrm < eps)
        L = 0;
        return
    end
    v = v / nrm;
end
% one more matvec to get the eigenvalue estimate
w = afun(v, 'notransp');
L = w(:).' * w(:);
end

function u = tv_prox_iso3d(z, lambda, niter)
% Chambolle's projected gradient for the isotropic 3D TV proximal:
%   argmin_u   0.5 * ||u - z||_2^2  +  lambda * TV_iso(u)
% Dual variables p = (p1, p2, p3); after the iterations,
%   u = z - lambda * div(p).

if (lambda <= 0)   % no TV penalty: the proximal operator is the identity
    u = z;
    return
end

sz = size(z);
p1 = zeros(sz);
p2 = zeros(sz);
p3 = zeros(sz);

tau = 1.0 / 12.0;       % ||D||_2^2 <= 12 in 3D forward differences

for k = 1:niter
    div_p = div3d(p1, p2, p3);
    aux = div_p - z / lambda;
    [g1, g2, g3] = grad3d(aux);

    p1 = p1 + tau * g1;
    p2 = p2 + tau * g2;
    p3 = p3 + tau * g3;

    % project onto the isotropic unit ball ||p_i|| <= 1 per voxel
    nrm = sqrt(p1 .* p1 + p2 .* p2 + p3 .* p3);
    nrm = max(nrm, 1.0);
    p1 = p1 ./ nrm;
    p2 = p2 ./ nrm;
    p3 = p3 ./ nrm;
end

u = z - lambda * div3d(p1, p2, p3);
end

function v = tv_iso3d_value(x)
% Isotropic TV value (used for objective history).
[g1, g2, g3] = grad3d(x);
v = sqrt(g1 .* g1 + g2 .* g2 + g3 .* g3);
v = sum(v(:));   % sum(..., 'all') is not supported by Octave < 7
end

function [g1, g2, g3] = grad3d(x)
% Forward differences with zero-Neumann BC (zero at last voxel along
% each axis).  Returns three arrays of the same shape as x.
g1 = zeros(size(x));
g2 = zeros(size(x));
g3 = zeros(size(x));
g1(1:end - 1, :, :) = x(2:end, :, :) - x(1:end - 1, :, :);
g2(:, 1:end - 1, :) = x(:, 2:end, :) - x(:, 1:end - 1, :);
g3(:, :, 1:end - 1) = x(:, :, 2:end) - x(:, :, 1:end - 1);
end

function d = div3d(p1, p2, p3)
% Discrete divergence (negative adjoint of grad3d).  Boundary handling
% matches grad3d so that <grad(x), p> = -<x, div(p)> exactly on the
% interior; the zero-Neumann BC pushes the boundary residue to zero.
d1 = zeros(size(p1));
d2 = zeros(size(p2));
d3 = zeros(size(p3));

d1(2:end - 1, :, :) = p1(2:end - 1, :, :) - p1(1:end - 2, :, :);
d1(1, :, :) = p1(1, :, :);
d1(end, :, :) = -p1(end - 1, :, :);

d2(:, 2:end - 1, :) = p2(:, 2:end - 1, :) - p2(:, 1:end - 2, :);
d2(:, 1, :) = p2(:, 1, :);
d2(:, end, :) = -p2(:, end - 1, :);

d3(:, :, 2:end - 1) = p3(:, :, 2:end - 1) - p3(:, :, 1:end - 2);
d3(:, :, 1) = p3(:, :, 1);
d3(:, :, end) = -p3(:, :, end - 1);

d = d1 + d2 + d3;
end
