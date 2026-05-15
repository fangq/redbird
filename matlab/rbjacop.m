function afun = rbjacop(J)
%
% afun = rbjacop(J)
%
% Wrap a 4D voxel-grid Jacobian as a function handle compatible with
% MATLAB's lsqr/lsmr:
%
%     afun(x, 'notransp') = J  * x        (matvec: voxel -> data)
%     afun(x, 'transp')   = J' * x        (rmatvec: data -> voxel)
%
% The 4D Jacobian shape is (Nx, Ny, Nz, Nsd), where Nsd = Ns*Nd is the
% number of source-detector pairs.  J(:,:,:,k) is the sensitivity volume
% for the k-th pair.  We never form the Nv x Nv Hessian J'*J; LSQR only
% needs the matvec/rmatvec products.
%
% For a 100^3 grid with 64 source-detector pairs:
%     J        :   ~250 MB   (the operator simply views this array)
%     J'*J     :   ~250 TB   (would never fit -- this is why LSQR is used)
%     LSQR cost:   2 mat-vec products per iteration, ~50-200 iters total.
%
% Usage:
%     afun = rbjacop(J);
%     delta = lsqr(afun, r, tol, maxit);   % r = y_meas - y_model
%
% Pair this with rbreglsqr for the full damped/early-stopped solve.
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     J: 4D Jacobian (Nx, Ny, Nz, Nsd) -- the voxel-grid adjoint
%        Jacobian returned by mcxlab in cfg.outputtype='adjoint'/
%        'adjoint_mua_d' mode.
%
% output:
%     afun: function handle with signature afun(x, flag) that returns:
%             - J * x        if flag = 'notransp'  (x: Nv x 1 -> Nsd x 1)
%             - J' * x       if flag = 'transp'    (x: Nsd x 1 -> Nv x 1)
%
% license:
%     GPL version 3, see LICENSE_GPLv3.txt files for details
%
% -- this function is part of Redbird-m toolbox
%

if (ndims(J) ~= 4)
    error('rbjacop: J must be a 4D array (Nx, Ny, Nz, Nsd)');
end

sz = size(J);
Nv = prod(sz(1:3));
Nsd = sz(4);

% Reshape J to a 2D (Nv x Nsd) view ONCE; both matvec and rmatvec then
% reduce to BLAS GEMV.  No copy made -- reshape on a contiguous array is
% an O(1) view in MATLAB.
J2 = reshape(J, Nv, Nsd);

afun = @(x, flag) apply(J2, Nv, Nsd, x, flag);

end

function y = apply(J2, Nv, Nsd, x, flag)
if (strcmp(flag, 'notransp'))
    if (numel(x) ~= Nv)
        error('rbjacop matvec: input length %d != Nv = %d', numel(x), Nv);
    end
    y = J2.' * x(:);
else
    if (numel(x) ~= Nsd)
        error('rbjacop rmatvec: input length %d != Nsd = %d', numel(x), Nsd);
    end
    y = J2 * x(:);
end
end
