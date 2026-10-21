function J_c = rbbinjac(J, binsize)
%
% J_c = rbbinjac(J, binsize)
%
% Bin a 4D fine-grid voxel Jacobian into a coarse-grid Jacobian by
% SUMMING fine-voxel columns within each coarse block.
%
% Mathematical interpretation: if x_c is the coarse-grid unknown and
% U is the upsample-by-replication operator that maps the coarse field
% to the fine field
%     x_fine[fine_voxel] = x_coarse[coarse_block(fine_voxel)],
% then the binned Jacobian is
%     J_c = J_fine * U
% i.e. each coarse Jacobian column = SUM of the binsize^3 fine
% Jacobian columns whose fine voxels lie in that coarse block.
%
% Why this is useful: voxel-grid MC adjoint Jacobians have one column
% per fine voxel (e.g. Nv = 1.08e5 for a 60 x 60 x 30 grid), making
% the system grossly underdetermined when Nsd <= a few hundred and
% forcing matrix-free iterative solvers (LSQR / FISTA-TV).  Binning to
% binsize = 2..6 yields Nv_c = a few hundred to a few thousand
% unknowns, comparable to a FEM mesh's node count, so the GN normal
% equations (J_c' * J_c + lambda * I) can be formed and solved
% directly -- exactly as in the original mesh-based redbird recon.
% Combine with rbupsample to map the coarse update back to the fine
% grid for the next forward call.
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     J:       4D fine Jacobian (Nx, Ny, Nz, Nsd), e.g. the .mua field
%              of the Jext struct returned by rbrunforward's mcxlab
%              adjoint path.
%     binsize: scalar or [bx, by, bz] integer bin factor.  Fine grids
%              whose dimensions are not exact multiples of binsize are
%              zero-padded on the trailing edge; round-trip with
%              rbupsample(..., fine_dims) trims the pad away.
%
% output:
%     J_c:     4D coarse Jacobian (Nx_c, Ny_c, Nz_c, Nsd) where
%              Nx_c = ceil(Nx / bx), etc.  Sd-axis is preserved.
%
% example:
%     J_c = rbbinjac(Jext.mua, 2);
%     [dmu_c, info] = rbreglsqr(J_c, r, 'maxit', 30, 'damp', 0.01);
%     dmu = rbupsample(dmu_c, 2, [Nx, Ny, Nz]);
%     mua_est = mua_est + single(dmu);
%
% see also: rbupsample, rbjacop, rbreglsqr, rbregtv
%
% license:
%     GPL version 3, see LICENSE_GPLv3.txt files for details
%
% -- this function is part of Redbird toolbox
%

if (ndims(J) ~= 4)
    error('rbbinjac: J must be a 4D array (Nx, Ny, Nz, Nsd)');
end
if (isscalar(binsize))
    binsize = [binsize, binsize, binsize];
end
if (any(binsize < 1) || any(binsize ~= round(binsize)))
    error('rbbinjac: binsize must be a positive integer or 3-vector');
end

sz = size(J);
Nx = sz(1);
Ny = sz(2);
Nz = sz(3);
Nsd = sz(4);
bx = binsize(1);
by = binsize(2);
bz = binsize(3);

% pad to next multiple of binsize so the reshape-and-sum trick lines up
Nx_p = bx * ceil(Nx / bx);
Ny_p = by * ceil(Ny / by);
Nz_p = bz * ceil(Nz / bz);

if (Nx_p ~= Nx || Ny_p ~= Ny || Nz_p ~= Nz)
    J_pad = zeros(Nx_p, Ny_p, Nz_p, Nsd, 'like', J);
    J_pad(1:Nx, 1:Ny, 1:Nz, :) = J;
else
    J_pad = J;
end

Nx_c = Nx_p / bx;
Ny_c = Ny_p / by;
Nz_c = Nz_p / bz;

% reshape into 7D so each bin axis is its own dimension, then collapse
% the three bin axes by summation.  MATLAB's column-major layout means
% the first bx samples of J_pad(:) span the first x-block, which is
% exactly what we want.
J_r = reshape(J_pad, bx, Nx_c, by, Ny_c, bz, Nz_c, Nsd);
J_r = sum(J_r, 1);
J_r = sum(J_r, 3);
J_r = sum(J_r, 5);

J_c = reshape(J_r, Nx_c, Ny_c, Nz_c, Nsd);
