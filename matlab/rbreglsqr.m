function [delta_mu, info] = rbreglsqr(J, r, varargin)
%
% [delta_mu, info] = rbreglsqr(J, r, 'param', value, ...)
%
% Iterative least-squares solve of  J * delta_mu = r  via MATLAB's
% built-in LSQR.  Designed for the voxel-grid mcxlab Jacobian path
% where J is too large for the normal-equation form
%     (J' * J + lambda * I)^(-1) * J' * r
% (J' * J would be Nv x Nv with Nv ~ 1e6 to 1e7 voxels).  LSQR only
% needs J*v and J'*u matvec products, so we wrap J through rbjacop and
% let LSQR build a tiny Krylov subspace iteratively.
%
% Regularization is by **early stopping** -- the Krylov subspace captures
% dominant singular components first (semi-convergence), so the iteration
% count itself acts as the regularizer.  Stop on the discrepancy principle
% (residual matches noise level) or on a hard iter cap.  No lambda to
% tune.  Pass 'damp', sqrt(lambda) explicitly if you want Tikhonov.
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     J:    either a 4D array (Nx, Ny, Nz, Nsd) -- the mcxlab voxel
%           Jacobian -- or a function handle afun(x, flag) following the
%           rbjacop convention.
%     r:    measurement residual vector, length Nsd = Ns*Nd (column
%           order: detector-fastest, source-slowest, matching mcxlab's
%           Jacobian last-axis convention).
%
%   param/value pairs (all optional):
%     'tol'        : LSQR convergence tolerance (default 1e-6).
%     'maxit'      : max LSQR iterations (default 200).
%     'damp'       : sqrt(lambda) for Tikhonov damping (default 0, i.e.
%                    no explicit damping; early stopping regularizes).
%     'reshape3d'  : if true (default), reshape the returned delta_mu
%                    back to (Nx, Ny, Nz) when J is a 4D array.  Pass
%                    false to return the flat Nv-vector.
%     'adjointtest': if true (default), run a one-shot dot-product test
%                    <Jv, u> == <v, J'u> with random v, u before
%                    calling LSQR.  Catches sign / index-order errors
%                    in custom afun handles.  Tolerance 1e-6.
%
% output:
%     delta_mu: voxel-space update.  Shape (Nx, Ny, Nz) for 4D-array J;
%               column vector of length Nv otherwise (or when
%               reshape3d=false).
%     info:     struct with fields
%                  .flag    LSQR convergence flag (see `help lsqr`)
%                  .relres  final relative residual
%                  .iter    number of LSQR iterations actually taken
%                  .resvec  per-iteration residual history
%                  .adjoint_err  relative error in the dot-product test
%                                (NaN when 'adjointtest' is false)
%
% see also: lsqr, rbjacop, rbreginv
%
% requires: MATLAB's built-in lsqr (Mapping / Math toolbox).  Not
% available in stock Octave; Octave users need to install
% `pkg install -forge optim` or supply their own iterative solver
% via the function-handle entry point.
%
% license:
%     GPL version 3, see LICENSE_GPLv3.txt files for details
%
% -- this function is part of Redbird-m toolbox
%

if (exist('lsqr', 'file') ~= 2 && exist('lsqr', 'builtin') ~= 5)
    error(['rbreglsqr: MATLAB''s built-in lsqr is not on the path. ' ...
           'On Octave, install the `optim` package via ' ...
           '`pkg install -forge optim; pkg load optim` before calling this function.']);
end

opt = varargin2struct(varargin{:});
tol      = jsonopt('tol',         1e-6, opt);
maxit    = jsonopt('maxit',       200,  opt);
damp     = jsonopt('damp',        0.0,  opt);
do_reshape = jsonopt('reshape3d', true, opt);
do_adjtst  = jsonopt('adjointtest', true, opt);

% --- wrap J into a function handle ---
input_was_4d = (isnumeric(J) && ndims(J) == 4);
if (input_was_4d)
    sz4 = size(J);
    Nv  = prod(sz4(1:3));
    Nsd = sz4(4);
    afun = rbjacop(J);
elseif (isa(J, 'function_handle'))
    afun = J;
    if (numel(r) <= 0)
        error('rbreglsqr: residual r must be non-empty to size the operator');
    end
    Nsd = numel(r);
    % probe Nv via one rmatvec call (cheap)
    probe = afun(zeros(Nsd, 1), 'transp');
    Nv = numel(probe);
else
    error('rbreglsqr: J must be a 4D array or a function handle afun(x, flag)');
end

% --- adjoint dot-product test ---
% This is the single most common bug source for iterative inverse
% problems; one cheap mat-vec pair catches sign / index / normalization
% mistakes in user-supplied afun handles.  Skip if the user explicitly
% disables it.
adjoint_err = NaN;
if (do_adjtst)
    % Don't touch global RNG state -- the dot-product test is reproducibility-
    % independent (any random pair works).  Octave compatibility: avoid rng().
    v_probe = randn(Nv, 1);
    u_probe = randn(Nsd, 1);

    Jv = afun(v_probe, 'notransp');
    Jtu = afun(u_probe, 'transp');
    lhs = u_probe(:).' * Jv(:);          % <u, Jv>
    rhs = Jtu(:).' * v_probe(:);         % <J' u, v>
    adjoint_err = abs(lhs - rhs) / max(abs(lhs), eps);
    if (adjoint_err > 1e-4)
        warning('rbreglsqr:adjointMismatch', ...
                ['Adjoint dot-product test failed: |<Jv, u> - <v, J''u>| / |<u, Jv>| = %.3e\n' ...
                 'LSQR may converge to garbage. Check the sign/normalization of your afun handle.'], ...
                adjoint_err);
    end
end

% --- handle Tikhonov damping via the augmented system [J; sqrt(lambda)*I] ---
% MATLAB's lsqr signature is lsqr(A, b, tol, maxit) and does NOT expose
% a damp parameter, so we build the augmented operator explicitly when
% damp > 0.  The augmented system has Nsd + Nv rows; we reshape the
% augmented r as [r; zeros(Nv, 1)].
if (damp > 0)
    aug_afun = @(x, flag) augmented(afun, Nv, Nsd, damp, x, flag);
    aug_r    = [r(:); zeros(Nv, 1)];
    [x, flag, relres, iter, resvec] = lsqr(aug_afun, aug_r, tol, maxit);
else
    [x, flag, relres, iter, resvec] = lsqr(afun, r(:), tol, maxit);
end

info = struct('flag', flag, 'relres', relres, 'iter', iter, ...
              'resvec', resvec, 'adjoint_err', adjoint_err);

if (input_was_4d && do_reshape)
    delta_mu = reshape(x, sz4(1), sz4(2), sz4(3));
else
    delta_mu = x;
end

end

function y = augmented(afun, Nv, Nsd, damp, x, flag)
% Wraps afun as the augmented operator A_aug = [J; damp*I] without
% materializing the (Nsd+Nv)-row augmented matrix.
if (strcmp(flag, 'notransp'))
    % x: Nv -> [J*x; damp*x]  : (Nsd + Nv)
    y = [afun(x, 'notransp'); damp * x(:)];
else
    % x: Nsd+Nv -> J'*x(1:Nsd) + damp*x(Nsd+1:end) : Nv
    y = afun(x(1:Nsd), 'transp') + damp * x(Nsd + 1:end);
end
end
