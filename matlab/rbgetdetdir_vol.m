function detdir = rbgetdetdir_vol(cfg)
%
% detdir=rbgetdetdir_vol(cfg)
%
% Estimate inward surface-normal directions at each detector position on
% a voxel-grid (cfg.vol) domain.
%
% This is the voxel-grid analogue of rbgetdetdir.m (which handles tet
% meshes).  The Monte Carlo (mcxlab) path in rbrunforward needs a
% Nd-by-4 cfg.detdir array (inward normal + focal length) to build
% detector-as-adjoint disk sources for the adjoint Jacobian computation.
% rbrunforward auto-invokes this helper when cfg.vol is given,
% cfg.detdir is missing, and a Jacobian is requested (nargout > 5).
%
% For each detector, the inward direction is estimated from the local
% gradient of the binary indicator (cfg.vol > 0) sampled in a small
% neighborhood.  This works for any domain shape whose boundary is
% reasonably resolved by the voxel grid; for very thin shells, increase
% the neighborhood half-width via cfg.detdirhw (default 3 voxels).
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     cfg: redbird simulation data structure with at minimum
%          cfg.vol  (Nx x Ny x Nz, non-zero inside the domain)
%          cfg.detpos (Nd x 3 or Nd x 4)
%          cfg.unitinmm (optional; defaults to 1)
%          cfg.detdirhw (optional; integer neighborhood half-width;
%                       default 3 voxels)
%
% output:
%     detdir: Nd-by-4 array.  Columns 1-3 are unit inward-normal
%             vectors; column 4 is the focal-length placeholder (0).
%
% license:
%     GPL version 3, see LICENSE_GPLv3.txt files for details
%
% -- this function is part of Redbird toolbox
%

if (~isfield(cfg, 'vol') || ~isfield(cfg, 'detpos'))
    error('rbgetdetdir_vol requires cfg.vol and cfg.detpos');
end

vol = cfg.vol;

% mcxlab's per-voxel-property modes (cfg.mediabyte > 100, e.g.
% MEDIA_MUA_FLOAT) pass cfg.vol as a 4D array (Nprop, Nx, Ny, Nz) where
% the first axis is the property channel.  For the indicator-based
% inward-normal estimate below we only need a 3D mask of inside-the-
% domain voxels; collapse the property axis.
if (ndims(vol) == 4)
    vol = squeeze(any(vol > 0, 1));
end
detpos = cfg.detpos;
unitinmm = 1.0;
if (isfield(cfg, 'unitinmm') && ~isempty(cfg.unitinmm))
    unitinmm = cfg.unitinmm;
end
hw = 3;
if (isfield(cfg, 'detdirhw') && ~isempty(cfg.detdirhw))
    hw = max(1, round(cfg.detdirhw));
end

[Nx, Ny, Nz] = size(vol);
inside = (vol > 0);

% Convert detector positions (in mm) into voxel-index coordinates.
% cfg.detpos uses the same 0-/1-based convention as mcxlab itself; here we
% treat it as 1-based grid coordinates -- the inward normal is rotation-
% independent so a one-voxel offset doesn't affect the direction.
detnum = size(detpos, 1);
detdir = zeros(detnum, 4);

for d = 1:detnum
    cx = round(detpos(d, 1) / unitinmm);
    cy = round(detpos(d, 2) / unitinmm);
    cz = round(detpos(d, 3) / unitinmm);

    % clamp to grid bounds
    cx = min(max(cx, 1), Nx);
    cy = min(max(cy, 1), Ny);
    cz = min(max(cz, 1), Nz);

    x0 = max(cx - hw, 1);
    x1 = min(cx + hw, Nx);
    y0 = max(cy - hw, 1);
    y1 = min(cy + hw, Ny);
    z0 = max(cz - hw, 1);
    z1 = min(cz + hw, Nz);

    block = double(inside(x0:x1, y0:y1, z0:z1));

    % vector field pointing FROM the detector voxel TO each block voxel,
    % weighted by the indicator: inside-voxels pull, outside-voxels push.
    [bx, by, bz] = ndgrid(x0:x1, y0:y1, z0:z1);
    dx = bx - cx;
    dy = by - cy;
    dz = bz - cz;
    rmag = sqrt(dx.^2 + dy.^2 + dz.^2);
    rmag(rmag == 0) = 1;   % skip self
    w = (2 * block - 1) ./ rmag;   % +1 for inside, -1 for outside

    nrm_x = sum(dx(:) .* w(:));
    nrm_y = sum(dy(:) .* w(:));
    nrm_z = sum(dz(:) .* w(:));

    nlen = sqrt(nrm_x^2 + nrm_y^2 + nrm_z^2);
    if (nlen < 1e-12)
        % degenerate (detector deep inside the volume); fall back to
        % the centroid heuristic used by rbgetdetdir on meshes
        domain_ctr = [(Nx + 1) / 2, (Ny + 1) / 2, (Nz + 1) / 2];
        v = domain_ctr - [cx, cy, cz];
        nlen = max(norm(v), 1e-12);
        detdir(d, 1:3) = v / nlen;
    else
        detdir(d, 1:3) = [nrm_x, nrm_y, nrm_z] / nlen;
    end
end
