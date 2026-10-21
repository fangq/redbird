function x_fine = rbupsample(x_coarse, binsize, fine_dims)
%
% x_fine = rbupsample(x_coarse, binsize, fine_dims)
%
% Upsample a 3D coarse field to a fine grid by REPLICATION.  Each
% coarse voxel value is copied to all binsize^3 fine voxels that lie
% in its block.
%
% This is the linear adjoint of the coarse-binning operation
% implemented by rbbinjac: with J_c = J_fine * U where U is the
% replication-upsample, rbupsample applies U to a coarse update so the
% resulting fine-grid delta_mu can be added back into the per-voxel
% mcxlab cfg.vol (or any per-voxel quantity).
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     x_coarse:  3D coarse-grid field, shape (Nx_c, Ny_c, Nz_c).
%     binsize:   scalar or [bx, by, bz] (must match the bin factor
%                that was used in rbbinjac).
%     fine_dims: optional [Nx, Ny, Nz] target fine-grid shape.  When
%                the original fine grid wasn't an exact multiple of
%                binsize, rbbinjac zero-padded the trailing edge --
%                pass fine_dims here to trim the upsample result back
%                to the original size.  Default: no trimming, output
%                shape is [Nx_c*bx, Ny_c*by, Nz_c*bz].
%
% output:
%     x_fine:    3D fine-grid field with each coarse value replicated
%                across its binsize^3 fine voxels.
%
% see also: rbbinjac, rbreglsqr, rbregtv
%
% license:
%     GPL version 3, see LICENSE_GPLv3.txt files for details
%
% -- this function is part of Redbird toolbox
%

if (isscalar(binsize))
    binsize = [binsize, binsize, binsize];
end

x_fine = repelem(x_coarse, binsize(1), binsize(2), binsize(3));

if (nargin >= 3 && ~isempty(fine_dims))
    sz_now = size(x_fine);
    if (any(sz_now(1:3) < fine_dims(1:3)))
        error('rbupsample: upsampled shape is smaller than requested fine_dims');
    end
    x_fine = x_fine(1:fine_dims(1), 1:fine_dims(2), 1:fine_dims(3));
end
