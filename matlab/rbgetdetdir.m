function detdir = rbgetdetdir(cfg)
%
% detdir=rbgetdetdir(cfg)
%
% Estimate inward surface-normal directions at each detector position.
%
% This helper exists for the Monte Carlo (mmclab) path in rbrunforward,
% which requires a Nd-by-4 cfg.detdir array (inward normal + focal length)
% to build detector-as-adjoint sources. rbrunforward auto-invokes
% rbgetdetdir when cfg.nphoton is set, cfg.detdir is missing, and a
% Jacobian is requested (nargout > 5 from rbrunforward).
%
% For each detector position, the nearest surface-face center is located
% and the face normal is oriented toward the mesh centroid (negated to
% produce the inward direction). This works well for convex or
% quasi-convex domains; for highly non-convex shapes (e.g. branching
% tubes) the user should provide cfg.detdir explicitly.
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     cfg: redbird simulation data structure with at minimum
%          cfg.node (Nn x 3), cfg.elem (Ne x 4 or 5), cfg.detpos (Nd x 3
%          or Nd x 4). If cfg.face is absent, it is built via volface.
%
% output:
%     detdir: Nd-by-4 array. Columns 1-3 are unit inward-normal vectors;
%             column 4 is the focal-length placeholder (0).
%
% license:
%     GPL version 3, see LICENSE_GPLv3.txt files for details
%
% -- this function is part of Redbird-m toolbox
%

if (~isfield(cfg, 'node') || ~isfield(cfg, 'elem') || ~isfield(cfg, 'detpos'))
    error('rbgetdetdir requires cfg.node, cfg.elem, and cfg.detpos');
end

if (~isfield(cfg, 'face') || isempty(cfg.face))
    cfg.face = volface(cfg.elem(:, 1:4));
end

% face normals - sign here depends on cfg.face vertex ordering; we orient
% below using the mesh centroid heuristic
ab = cfg.node(cfg.face(:, 2), 1:3) - cfg.node(cfg.face(:, 1), 1:3);
ac = cfg.node(cfg.face(:, 3), 1:3) - cfg.node(cfg.face(:, 1), 1:3);
nrm = cross(ab, ac, 2);
nlen = sqrt(sum(nrm .* nrm, 2));
facenrm = nrm ./ repmat(nlen, 1, 3);

% face center coordinates
facectr = (cfg.node(cfg.face(:, 1), 1:3) + ...
           cfg.node(cfg.face(:, 2), 1:3) + ...
           cfg.node(cfg.face(:, 3), 1:3)) / 3;

% orient each face normal to point OUTWARD (away from mesh centroid),
% then negate to get the inward direction used by mmc as the adjoint
% source emission direction
meshctr = mean(cfg.node(:, 1:3), 1);
out_proj = sum(facenrm .* (facectr - repmat(meshctr, size(facectr, 1), 1)), 2);
flipsign = sign(out_proj);
flipsign(flipsign == 0) = 1;
facenrm = facenrm .* repmat(flipsign, 1, 3);

% nearest-face lookup for each detector
detnum = size(cfg.detpos, 1);
detdir = zeros(detnum, 4);
for d = 1:detnum
    delta = facectr - repmat(cfg.detpos(d, 1:3), size(facectr, 1), 1);
    d2 = sum(delta .* delta, 2);
    [~, fid] = min(d2);
    detdir(d, 1:3) = -facenrm(fid, :);
end
