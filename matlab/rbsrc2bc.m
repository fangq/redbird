function cfg = rbsrc2bc(cfg, isdet)
%
% srcbc=rbsrc2bc(cfg)
%
% Converting wide-field source forms into a boundary condition by defining
% in-ward flux on the mesh surface
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     cfg: the simulation data structure, with srctype, srcpos, srcdir,
%          srcparam1, srcparam2, srcpattern fields
%     isdet: default is 0; if set to 1, rbsrc2bc process widefield
%          detectors, the relevant fields are dettype, detpos, detdir,
%          detparam1, detparam2, detpattern
%
% output:
%     srcbc: an array of Ns x Nt, where Nt is size(cfg.face,1) and Ns is
%            the number of sources (if isdet=1, the detector counts)
%
% license:
%     GPL version 3, see LICENSE_GPLv3.txt files for details
%
% -- this function is part of Redbird-m toolbox
%

srcbc = [];

if (nargin < 2)
    isdet = 0;
end

% line-source path (works for both DOT and MWT): triggered by
% cfg.srctype=='line' (or cfg.dettype=='line' for detectors). cfg.srcpos(:,1:3)
% is the line start; cfg.srcpos(:,1:3)+cfg.srcparam1(1:3) is the end. Each
% line is traced through the tet mesh and contributes a sparse RHS column
% built from per-tet line integrals of the linear basis. No 1/musp sinking
% is applied (unlike pencil sources).
if (~isdet && isfield(cfg, 'srctype') && strcmp(cfg.srctype, 'line') && isfield(cfg, 'srcpos') && ~isempty(cfg.srcpos))
    if (~isfield(cfg, 'srcparam1') || numel(cfg.srcparam1) < 3)
        error('cfg.srctype=''line'' requires cfg.srcparam1(1:3) (line endpoint offset)');
    end
    cfg.srcpos0 = cfg.srcpos;
    linesegs = [cfg.srcpos(:, 1:3), cfg.srcpos(:, 1:3) + repmat(cfg.srcparam1(1:3), size(cfg.srcpos, 1), 1)];
    cfg.widesrc = lineseg2rhs(linesegs, cfg);
    cfg.srcpos = zeros(0, 3);
    return
end
if (isdet && isfield(cfg, 'dettype') && strcmp(cfg.dettype, 'line') && isfield(cfg, 'detpos') && ~isempty(cfg.detpos))
    if (~isfield(cfg, 'detparam1') || numel(cfg.detparam1) < 3)
        error('cfg.dettype=''line'' requires cfg.detparam1(1:3) (line endpoint offset)');
    end
    cfg.detpos0 = cfg.detpos;
    linesegs = [cfg.detpos(:, 1:3), cfg.detpos(:, 1:3) + repmat(cfg.detparam1(1:3), size(cfg.detpos, 1), 1)];
    cfg.widedet = lineseg2rhs(linesegs, cfg);
    cfg.detpos = zeros(0, 3);
    return
end

% ray-source path: cfg.srctype=='ray' models a collimated beam as a line of
% isotropic sources with weight density mu_s' * exp(-mu_s' * l), normalized so
% the total integral over the (untruncated) ray equals 1 and the mean depth is
% 1/mu_s' = l_tr (Haskell et al. 1994, eq 2.4.5). Truncated at amplitude
% threshold (default 1e-3) along cfg.srcdir from cfg.srcpos. cfg.srcdir is
% required (unlike 'line' which uses srcparam1). Same applies for 'ray' dets.
if (~isdet && isfield(cfg, 'srctype') && strcmp(cfg.srctype, 'ray') && isfield(cfg, 'srcpos') && ~isempty(cfg.srcpos))
    if (~isfield(cfg, 'srcdir') || isempty(cfg.srcdir))
        error('cfg.srctype=''ray'' requires cfg.srcdir');
    end
    cfg.srcpos0 = cfg.srcpos;
    cfg.widesrc = rayseg2rhs(cfg.srcpos, cfg.srcdir, cfg);
    cfg.srcpos = zeros(0, 3);
    return
end
if (isdet && isfield(cfg, 'dettype') && strcmp(cfg.dettype, 'ray') && isfield(cfg, 'detpos') && ~isempty(cfg.detpos))
    if (~isfield(cfg, 'detdir') || isempty(cfg.detdir))
        error('cfg.dettype=''ray'' requires cfg.detdir');
    end
    cfg.detpos0 = cfg.detpos;
    cfg.widedet = rayseg2rhs(cfg.detpos, cfg.detdir, cfg);
    cfg.detpos = zeros(0, 3);
    return
end

if (~isdet)
    if ((~isfield(cfg, 'srctype') && ~isfield(cfg, 'widesrcid')) || (isfield(cfg, 'srctype') && (strcmp(cfg.srctype, 'pencil') || strcmp(cfg.srctype, 'isotropic'))))
        return
    end
    srcdir = cfg.srcdir;
    sources = cfg.srcpos;
    if isfield(cfg, 'widesrcid')
        widesrcid = cfg.widesrcid;
        if ~isa(widesrcid, 'containers.Map')
            widesrcid = containers.Map({'_'}, {widesrcid});
        end
    else
        tempwf.srctype = {cfg.srctype};
        if isfield(cfg, 'srcid')
            tempwf.srcid = {cfg.srcid};
        else
            tempwf.srcid = {1};
        end
        try
            tempwf.srcparam1 = {cfg.srcparam1};
            tempwf.srcparam2 = {cfg.srcparam2};
        catch
            error('Widefield sources must have srcparam1 and srcparam2');
        end
        if isfield(cfg, 'srcpattern')
            tempwf.srcpattern = {cfg.srcpattern};
        end
        if isfield(cfg, 'srcweight')
            tempwf.srcweight = {cfg.srcweight};
        end
        widesrcid = containers.Map({'_'}, {tempwf});
        clear tempwf;
    end
else
    if ((~isfield(cfg, 'dettype') && ~isfield(cfg, 'widedetid')) || (isfield(cfg, 'dettype') && (strcmp(cfg.dettype, 'pencil') || strcmp(cfg.dettype, 'isotropic'))))
        return
    end
    srcdir = cfg.detdir;
    sources = cfg.detpos;
    if isfield(cfg, 'widedetid')
        widesrcid = cfg.widedetid;
        if ~isa(widesrcid, 'containers.Map')
            widesrcid = containers.Map({'_'}, {widesrcid});
        end
    else
        tempwf.srctype = {cfg.dettype};
        if isfield(cfg, 'detid')
            tempwf.srcid = {cfg.detid};
        else
            tempwf.srcid = {1};
        end
        try
            tempwf.srcparam1 = {cfg.detparam1};
            tempwf.srcparam2 = {cfg.detparam2};
        catch
            error('Widefield detectors must have detparam1 and detparam2');
        end
        if isfield(cfg, 'detpattern')
            tempwf.srcpattern = {cfg.detpattern};
        end
        if isfield(cfg, 'detweight')
            tempwf.srcweight = {cfg.detweight};
        end
        widesrcid = containers.Map({'_'}, {tempwf});
        clear tempwf;
    end
end

if isa(cfg.prop, 'containers.Map')
    prop = cfg.prop;
else
    prop = containers.Map({'_'}, {cfg.prop});
end

if (~isequal(cell2mat(widesrcid.keys), cell2mat(prop.keys)) && (~isempty(cell2mat(widesrcid.keys))))
    error('Keys for prop and widesrcid must be compatible');
elseif (~isequal(cell2mat(widesrcid.keys), cell2mat(prop.keys)) && (~isempty(cell2mat(prop.keys))))
    temp = widesrcid('');
    clear widesrcid;
    widesrcid = containers.Map();
    wv = prop.keys;
    for ii = 1:length(prop.keys)
        widesrcid(wv{ii}) = temp;
    end
end

widesrc = [];
wavelengths = prop.keys;
wfsrcmapping = containers.Map();
allidwf = [];

for wv = wavelengths
    wideparam = widesrcid(wv{1});
    op = prop(wv{1});

    z0 = 1 / (op(2, 1) + op(2, 2) * (1 - op(2, 3)));

    srcmapping = [];

    for wideidx = 1:length(wideparam.srcid)
        clear srcbc rhs;

        srcid = wideparam.srcid{wideidx};
        allidwf = [allidwf; srcid];
        srctype = wideparam.srctype{wideidx};
        srcparam1 = wideparam.srcparam1{wideidx};
        srcparam2 = wideparam.srcparam2{wideidx};

        if (strcmp(srctype, 'pattern'))
            srcpattern = wideparam.srcpattern{wideidx};
        end
        if isfield(wideparam, 'srcweight')
            srcweight = wideparam.srcweight{wideidx};
        end

        srcpos = sources(srcid, :);

        switch srctype
            case {'planar', 'pattern', 'fourier'}
                ps = [srcpos; srcpos + srcparam1(1:3); ...
                      srcpos + srcparam1(1:3) + srcparam2(1:3); srcpos + srcparam2(1:3); srcpos];

                pnode = cfg.node;
                pface = cfg.face;

                % if src is colimated (default), sink it by 1/mus'
                if (~isfield(cfg, 'iscolimated') || cfg.iscolimated)
                    sinkplane = srcdir; % the plane where sinked planar source is located as [A,B,C,D] where A*x+B*y+C*z+D=0
                    sinkplane(4) = -sum(srcdir .* (srcpos + srcdir * z0));
                    [cutpos, cutvalue, facedata, elemid, nodeid] = qmeshcut(cfg.elem, cfg.node, zeros(size(cfg.node, 1), 1), sinkplane);
                    pnode = cutpos;
                    idx = find(facedata(:, 3) ~= facedata(:, 4));
                    pface = facedata(facedata(:, 3) == facedata(:, 4), 1:3);
                    pface = [pface; facedata(idx, [1 2 3]); facedata(idx, [1 3 4])];
                end
                c0 = meshcentroid(pnode, pface);
                newnode = rotatevec3d([c0; ps], srcdir(1:3));
                srcpoly = newnode(end - 4:end, 1:2);
                [isin, ison] = inpolygon(newnode(1:end - 5, 1), newnode(1:end - 5, 2), srcpoly(:, 1), srcpoly(:, 2));
                isin = isin | ison;
                idx = find(isin);
                if (~isempty(idx)) % the below test only works for convex shapes
                    AB = pnode(pface(idx, 2), 1:3) - pnode(pface(idx, 1), 1:3);
                    AC = pnode(pface(idx, 3), 1:3) - pnode(pface(idx, 1), 1:3);
                    N = cross(AB', AC')';
                    dir = sum(N .* repmat(srcdir(:)', size(N, 1), 1), 2);
                    if (exist('sinkplane', 'var'))
                        dir(dir > 0) = -dir(dir > 0);
                    end
                    if (all(dir >= 0))
                        error('please reorient the surface triangles');
                    end
                    srcbc = zeros(1, size(pface, 1));
                    srcbc(idx(dir < 0)) = 1;
                    pbc = newnode(idx(dir < 0), 1:2);

                    dp = pbc - repmat(srcpoly(1, :), size(pbc, 1), 1);
                    dx = srcpoly(2, :) - srcpoly(1, :);
                    dy = srcpoly(4, :) - srcpoly(1, :);
                    nx = dx / norm(dx);
                    ny = dy / norm(dy);

                    bary = [sum(dp .* repmat(nx / norm(dx), size(dp, 1), 1), 2), sum(dp .* repmat(ny / norm(dy), size(dp, 1), 1), 2)];
                    bary(bary < 0) = 0;
                    bary(bary >= 1) = 1 - 1e-6;

                    if (exist('srcpattern', 'var'))
                        srcpattern = permute(srcpattern, [3 1 2]);
                        pdim = size(srcpattern);
                        patsize = pdim(1);
                        srcbc = repmat(srcbc, patsize, 1);
                        for i = 1:patsize
                            srcbc(i, idx(dir < 0)) = srcpattern(sub2ind(pdim, i * ones(size(bary, 1), 1), floor(bary(:, 1) * pdim(2)) + 1, floor(bary(:, 2) * pdim(3)) + 1));
                        end
                    elseif (strcmp(srctype, 'fourier'))
                        kx = floor(srcparam1(4));
                        ky = floor(srcparam2(4));
                        phi0 = (srcparam1(4) - kx) * 2 * pi;
                        M = 1 - (srcparam2(4) - ky);
                        srcbc = repmat(srcbc, kx * ky, 1);
                        for i = 1:kx
                            for j = 1:ky
                                srcbc((i - 1) * ky + j, idx(dir < 0)) = 0.5 * (1 + M * cos((i * bary(:, 1) + j * bary(:, 2)) * 2 * pi + phi0))';
                            end
                        end
                    end
                else
                    error('source direction does not intersect with the domain');
                end
            otherwise
                error('this source type is not supported');
        end

        % at this point, srcbc stores the J- at each surface triangle (or sinked triangles)

        if isa(cfg.reff, 'containers.Map')
            Reff = cfg.reff(wv{1});
        else
            Reff = cfg.reff;
        end
        maxbcnode = max(pface(:));
        if (exist('nodeid', 'var'))
            nodeweight = nodeid(:, 3);
            nodeid = nodeid(:, 1:2);
            maxbcnode = max(max(nodeid(pface, :)));
            parea = elemvolume(pnode, pface);
        else
            parea = cfg.area;
        end

        % 1/18 = 1/2*1/9, where 2 comes from the 1/2 in ls=(1+Reff)/(1-Reff)/2*D,
        % and 1/9 = (1/6+1/12+1/12)/3, where A/6 is <phi_i,phi_j> when i=j, and
        % A/12 is i!=j
        Adiagbc = parea(:) * ((1 - Reff) / (18 * (1 + Reff))); % ones(size(parea(:),1),1)*((1-Reff)/(18*(1+Reff)));
        Adiagbc = repmat(Adiagbc, 1, size(srcbc, 1)) .* (srcbc');
        rhs = sparse(size(cfg.node, 1), size(srcbc, 1));

        for i = 1:size(srcbc, 1)
            if (exist('nodeid', 'var'))
                allnodes = nodeid(pface, :);
                rhs(1:maxbcnode, i) = sparse(allnodes(:), 1, [repmat(Adiagbc(:, i), 3, 1) .* nodeweight(pface(:)); repmat(Adiagbc(:, i), 3, 1) .* (1 - nodeweight(pface(:)))]);
            else
                rhs(1:maxbcnode, i) = sparse(cfg.face(:), 1, repmat(Adiagbc(:, i), 1, 3));
            end
            wsrc = 1;
            if (exist('srcweight', 'var') && numel(srcweight) == size(srcbc, 1))
                wsrc = srcweight(i);
            end
            %             rhs(rhs(:,i)<1e-3*max(rhs(:,i)),i) = 0;
            %             rhs(:,i)=rhs(:,i)*(wsrc/abs(sum(rhs(:,i))));
            rhs(:, i) = rhs(:, i) * (wsrc / sum(abs(rhs(:, i))));
            %             rhs(:,i)=rhs(:,i)*(wsrc/sum(rhs(:,i)));
        end
        srcbc = rhs.';
        if exist('srcpattern', 'var')
            patsize = size(srcpattern, 1);
        else
            patsize = 1;
        end

        indices = [(size(widesrc, 1) + 1) (size(widesrc, 1) + patsize)];
        widesrc = [widesrc; full(srcbc)];
        srcmapping = [srcmapping; srcid indices];
    end
    wfsrcmapping(wv{1}) = srcmapping;
end
sources(unique(allidwf), :) = [];

if (length(wfsrcmapping) == 1)
    wfsrcmapping = wfsrcmapping(wavelengths{1});
end

if (~isdet)
    cfg.widesrc = widesrc;
    cfg.wfsrcmapping = wfsrcmapping;
    cfg.srcpos = sources;
else
    cfg.widedet = widesrc;
    cfg.wfdetmapping = wfsrcmapping;
    cfg.detpos = sources;
end

function widesrc = lineseg2rhs(linesegs, cfg)
% build the RHS for each line segment source (columns of the result are the
% RHS contribution at all Nn nodes; one row per line). For Helmholtz, scale
% by -j*omega*mu0 (mu0 in H/mm).

ns = size(linesegs, 1);
nn = size(cfg.node, 1);
widesrc = zeros(ns, nn);

for s = 1:ns
    rowvec = traceline(linesegs(s, 1:3), linesegs(s, 4:6), cfg);
    widesrc(s, :) = rowvec;
end

if isfield(cfg, 'srcweight') && numel(cfg.srcweight) == ns
    widesrc = widesrc .* cfg.srcweight(:);
end

ishelmholtz = isfield(cfg, 'bulk') && (isfield(cfg.bulk, 'epsilon') || isfield(cfg.bulk, 'sigma'));
if (ishelmholtz)
    mu0_mm = 4 * pi * 1e-10;
    if isa(cfg.omega, 'containers.Map')
        kk = cfg.omega.keys;
        omega = cfg.omega(kk{1});
    else
        omega = cfg.omega;
    end
    widesrc = -1j * omega * mu0_mm * widesrc;
end

function rowvec = traceline(p1, p2, cfg)
% trace one line segment through the tet mesh and accumulate
% (bary_in + bary_out)/2 * dl per tet onto the per-node line integral.

nn = size(cfg.node, 1);
rowvec = zeros(1, nn);

p1 = p1(:)';
p2 = p2(:)';
seglen = norm(p2 - p1);
if seglen < 1e-12
    return
end
dirvec = (p2 - p1) / seglen;

[estart, barystart] = tsearchn(cfg.node(:, 1:3), cfg.elem(:, 1:4), p1);
if isnan(estart)
    return
end

face_local = [1 2 3; 1 2 4; 1 3 4; 2 3 4];

e = estart;
pcur = p1;
barycur = barystart(:)';
lrem = seglen;

maxiter = 100 * size(cfg.elem, 1);
iter = 0;
while lrem > 1e-12 && e > 0 && iter < maxiter
    iter = iter + 1;
    elemnodes = cfg.elem(e, 1:4);
    face_global = elemnodes(face_local);
    [tt, uu, vv] = raytrace(pcur, dirvec, cfg.node(:, 1:3), face_global);
    valid = (tt > 1e-10) & ~isinf(tt) & (uu >= -1e-10) & (vv >= -1e-10) & (uu + vv <= 1.0 + 1e-10);
    validid = find(valid);
    if isempty(validid)
        break
    end
    [tmin, kk] = min(tt(validid));
    flocal = validid(kk);

    if tmin >= lrem - 1e-12
        % segment ends inside this tet
        pend = pcur + lrem * dirvec;
        verts = cfg.node(elemnodes, 1:3)';
        baryend = ([verts; ones(1, 4)] \ [pend(:); 1])';
        contrib = (barycur + baryend) * 0.5 * lrem;
        rowvec(elemnodes) = rowvec(elemnodes) + contrib;
        break
    end

    u_e = uu(flocal);
    v_e = vv(flocal);
    fnodes = face_local(flocal, :);
    baryexit = zeros(1, 4);
    baryexit(fnodes(1)) = 1 - u_e - v_e;
    baryexit(fnodes(2)) = u_e;
    baryexit(fnodes(3)) = v_e;

    contrib = (barycur + baryexit) * 0.5 * tmin;
    rowvec(elemnodes) = rowvec(elemnodes) + contrib;

    enext = cfg.facenb(e, flocal);
    pcur = pcur + tmin * dirvec;
    lrem = lrem - tmin;
    e = enext;
    if e > 0
        % re-express bary in the new tet's local node order (the 3 shared face
        % nodes typically occupy different local positions in e vs enext)
        elemnodes_next = cfg.elem(e, 1:4);
        verts_next = cfg.node(elemnodes_next, 1:3)';
        barycur = ([verts_next; ones(1, 4)] \ [pcur(:); 1])';
    else
        barycur = baryexit;
    end
end

function widesrc = rayseg2rhs(positions, dirs, cfg)
% build the RHS for each ray source/detector. Each ray starts at positions(s,:)
% and propagates along dirs(s,:) (or dirs broadcast), with weight density
% mu_s' * exp(-mu_tr * l) along the beam (Haskell et al. 1994, eq 2.4.5).
% mu_tr = mu_a + mu_s' is the transport-attenuation coefficient (the photon
% is depleted by both scattering and absorption en route to the first-scatter
% point); mu_s' is the rate of scattering events per unit length. The total
% integral over the (untruncated) ray is the transport albedo mu_s'/mu_tr.
%
% Geometric tet-walk is performed once per ray (wavelength-independent); the
% per-tet records are then integrated using each wavelength's (mu_s', mu_tr)
% pair. For multi-wavelength runs (cfg.musp0 a containers.Map), the cached
% geometry can be reintegrated; this prototype outputs a single matrix using
% the first wavelength's coefficients. For Helmholtz/MWT, widesrc is scaled by
% -j*omega*mu0 to match the line-source convention.

ns = size(positions, 1);
nn = size(cfg.node, 1);

threshold = 1e-3;
if (isfield(cfg, 'raythreshold') && ~isempty(cfg.raythreshold))
    threshold = cfg.raythreshold;
end

% bulk optical properties to compute (mu_s', mu_tr) per wavelength
bkprop = rbgetbulk(cfg);
if isa(bkprop, 'containers.Map')
    keys = bkprop.keys;
    musp_vals  = zeros(1, length(keys));
    mutr_vals = zeros(1, length(keys));
    for k = 1:length(keys)
        bk = bkprop(keys{k});
        musp_vals(k)  = bk(2) * (1 - bk(3));
        mutr_vals(k) = bk(1) + musp_vals(k);
    end
    musp_default  = musp_vals(1);
    mutr_default = mutr_vals(1);
    mutr_min     = min(mutr_vals);
else
    musp_default  = bkprop(2) * (1 - bkprop(3));
    mutr_default = bkprop(1) + musp_default;
    mutr_min     = mutr_default;
end

if (isfield(cfg, 'raylen') && ~isempty(cfg.raylen))
    L_max = cfg.raylen;
else
    L_max = -log(threshold) / mutr_min;
end

if size(dirs, 1) == 1
    dirs = repmat(dirs, ns, 1);
end
dirnorms = sqrt(sum(dirs.^2, 2));
dirs = dirs ./ repmat(dirnorms, 1, 3);

widesrc = zeros(ns, nn);
for s = 1:ns
    geom = ray_trace_geom(positions(s, 1:3), dirs(s, :), L_max, cfg);
    widesrc(s, :) = ray_integrate(geom, musp_default, mutr_default, cfg);
end

if isfield(cfg, 'srcweight') && numel(cfg.srcweight) == ns
    widesrc = widesrc .* cfg.srcweight(:);
end

ishelmholtz = isfield(cfg, 'bulk') && (isfield(cfg.bulk, 'epsilon') || isfield(cfg.bulk, 'sigma'));
if (ishelmholtz)
    mu0_mm = 4 * pi * 1e-10;
    if isa(cfg.omega, 'containers.Map')
        kk = cfg.omega.keys;
        omega = cfg.omega(kk{1});
    else
        omega = cfg.omega;
    end
    widesrc = -1j * omega * mu0_mm * widesrc;
end

function geom = ray_trace_geom(p0, dirvec, L_max, cfg)
% walk the ray through tetrahedra up to length L_max; record per-tet info
% (element id, entry/exit distances along the ray, entry/exit barycentric
% coordinates). geometry is wavelength-independent and can be reused for
% multi-wavelength integration.

p0 = p0(:)';
dirvec = dirvec(:)';

% Tiny perpendicular nudge on p0 to break grid-aligned ray failures: when
% dirvec is axis-aligned and p0 lies on a face shared by 4+ structured-mesh
% tets, raytrace returns Inf for two faces and 0 for two faces of the
% post-transition tet, leaving no valid forward exit and silently aborting
% the trace after the 1st tet (mass loss ~30% for typical setups). The
% nudge magnitude is ~1e-6 mm, far below any mesh resolution, so physics is
% unchanged.
abs_dir = abs(dirvec);
[~, idx] = min(abs_dir);
e_arb = zeros(1, 3);
e_arb(idx) = 1;
perp1 = cross(dirvec, e_arb);
perp1 = perp1 / norm(perp1);
perp2 = cross(dirvec, perp1);
perp2 = perp2 / norm(perp2);
p0 = p0 + 1e-6 * perp1 + 1.3e-6 * perp2;

geom = struct('elemids', zeros(0, 1), 'zs', zeros(0, 2), ...
              'barys_in', zeros(0, 4), 'barys_out', zeros(0, 4));

[estart, barystart] = tsearchn(cfg.node(:, 1:3), cfg.elem(:, 1:4), p0);
if isnan(estart)
    return
end

face_local = [1 2 3; 1 2 4; 1 3 4; 2 3 4];

e = estart;
pcur = p0;
barycur = barystart(:)';
zcur = 0;
lrem = L_max;

elemids = [];
zs = [];
barys_in = [];
barys_out = [];

maxiter = 100 * size(cfg.elem, 1);
iter = 0;
while lrem > 1e-12 && e > 0 && iter < maxiter
    iter = iter + 1;
    elemnodes = cfg.elem(e, 1:4);
    face_global = elemnodes(face_local);
    [tt, uu, vv] = raytrace(pcur, dirvec, cfg.node(:, 1:3), face_global);
    valid = (tt > 1e-10) & ~isinf(tt) & (uu >= -1e-10) & (vv >= -1e-10) & (uu + vv <= 1.0 + 1e-10);
    validid = find(valid);
    if isempty(validid)
        break
    end
    [tmin, kk] = min(tt(validid));
    flocal = validid(kk);

    if tmin >= lrem - 1e-12
        % ray ends inside this tet (truncated by L_max)
        pend = pcur + lrem * dirvec;
        verts = cfg.node(elemnodes, 1:3)';
        baryend = ([verts; ones(1, 4)] \ [pend(:); 1])';
        elemids = [elemids; e];
        zs = [zs; zcur, zcur + lrem];
        barys_in = [barys_in; barycur];
        barys_out = [barys_out; baryend];
        break
    end

    u_e = uu(flocal);
    v_e = vv(flocal);
    fnodes = face_local(flocal, :);
    baryexit = zeros(1, 4);
    baryexit(fnodes(1)) = 1 - u_e - v_e;
    baryexit(fnodes(2)) = u_e;
    baryexit(fnodes(3)) = v_e;

    elemids = [elemids; e];
    zs = [zs; zcur, zcur + tmin];
    barys_in = [barys_in; barycur];
    barys_out = [barys_out; baryexit];

    enext = cfg.facenb(e, flocal);
    pcur = pcur + tmin * dirvec;
    zcur = zcur + tmin;
    lrem = lrem - tmin;
    e = enext;
    if e > 0
        % re-express bary in the new tet's local node order (face nodes are
        % shared but typically permuted in the new tet's local indexing)
        elemnodes_next = cfg.elem(e, 1:4);
        verts_next = cfg.node(elemnodes_next, 1:3)';
        barycur = ([verts_next; ones(1, 4)] \ [pcur(:); 1])';
    else
        barycur = baryexit;
    end
end

geom.elemids = elemids;
geom.zs = zs;
geom.barys_in = barys_in;
geom.barys_out = barys_out;

function rowvec = ray_integrate(geom, amp, rate, cfg)
% accumulate the closed-form integral amp * exp(-rate * z) * baryc(z) across
% all recorded tets onto an Nn-row vector. baryc is linear in z within each
% tet, so the per-tet integral is exact. amp and rate are decoupled so the
% physically motivated weights (amp = mu_s', rate = mu_tr) can be set
% independently. With amp = mu_s' and rate = mu_tr, total ray integral
% equals the transport albedo mu_s'/mu_tr (=1 only when mu_a -> 0).

nn = size(cfg.node, 1);
rowvec = zeros(1, nn);
ntet = size(geom.elemids, 1);
if ntet == 0
    return
end
for k = 1:ntet
    e = geom.elemids(k);
    z_in = geom.zs(k, 1);
    z_out = geom.zs(k, 2);
    b_in = geom.barys_in(k, :);
    b_out = geom.barys_out(k, :);
    dl = z_out - z_in;
    if dl <= 0
        continue
    end
    x = rate * dl;
    if abs(x) < 1e-3
        % Taylor expansion to avoid catastrophic cancellation when rate*dl ~ 0
        A = dl * (1 - x / 2 + x * x / 6 - x * x * x / 24);
        B = dl * dl * (0.5 - x / 3 + x * x / 8 - x * x * x / 30);
    else
        e_d = exp(-x);
        A = -expm1(-x) / rate;          % integral of exp(-rate*u) from 0 to dl
        B = A / rate - dl * e_d / rate; % integral of u*exp(-rate*u) from 0 to dl
    end
    scale = amp * exp(-rate * z_in);
    contrib = scale * (b_in * A + (b_out - b_in) * (B / dl));
    elemnodes = cfg.elem(e, 1:4);
    rowvec(elemnodes) = rowvec(elemnodes) + contrib;
end
