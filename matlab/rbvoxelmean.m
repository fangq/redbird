function v = rbvoxelmean(phi, pts, avgsize)
%
% v = rbvoxelmean(phi, pts, avgsize)
%
% Sample a volumetric field at a set of points by averaging voxel
% values within a cubic neighborhood of half-width `avgsize` (in
% voxels) around each point.
%
% Used by the rbrunforward mcxlab branch to extract boundary detector
% readings from a 3D fluence volume.  Nearest-neighbor + small box
% averaging is preferred over strict trilinear interpolation because
% detectors at the surface land on voxels whose neighbors are outside
% the domain (and therefore zero), which a trilinear scheme would
% halve.  Box averaging matches the conventional "detector aperture"
% interpretation used in DOT.
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     phi:     Nx x Ny x Nz volumetric field (e.g. fluence).
%     pts:     Np x 3 point coordinates, in voxel-index units (the
%              caller is expected to have applied unitinmm scaling
%              already, matching mcxlab's cfg.detpos convention).
%     avgsize: optional, scalar half-width in voxels of the averaging
%              box (default 1, i.e. 3 x 3 x 3 cube; pass 0 for pure
%              nearest-neighbor sampling).
%
% output:
%     v:      Np x 1 column vector of averaged values, one per point.
%
% license:
%     GPL version 3, see LICENSE_GPLv3.txt files for details
%
% -- this function is part of Redbird toolbox
%

if (nargin < 3 || isempty(avgsize))
    avgsize = 1;
end
hw = max(0, round(avgsize));

[Nx, Ny, Nz] = size(phi);
Np = size(pts, 1);
v = zeros(Np, 1);

for k = 1:Np
    cx = round(pts(k, 1));
    cy = round(pts(k, 2));
    cz = round(pts(k, 3));

    cx = min(max(cx, 1), Nx);
    cy = min(max(cy, 1), Ny);
    cz = min(max(cz, 1), Nz);

    if (hw == 0)
        v(k) = phi(cx, cy, cz);
        continue
    end

    x0 = max(cx - hw, 1);
    x1 = min(cx + hw, Nx);
    y0 = max(cy - hw, 1);
    y1 = min(cy + hw, Ny);
    z0 = max(cz - hw, 1);
    z1 = min(cz + hw, Nz);

    block = phi(x0:x1, y0:y1, z0:z1);
    v(k) = mean(block(:));
end
