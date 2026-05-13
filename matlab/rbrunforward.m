function [detval, phi, Amat, rhs, sflag, Jext] = rbrunforward(cfg, varargin)
%
% [detval, phi]=rbrunforward(cfg)
%    or
% [detval, phi, Amat, rhs]=rbrunforward(cfg,'param1',value1,...)
%    or
% [detval, phi, Amat, rhs, sflag, Jext]=rbrunforward(cfg,'param1',value1,...)
%
% Perform forward simulations at all sources and all wavelengths based on the input structure
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     cfg: the redbird forward simulation data structure
%       The forward solver in Redbird follows a very similar structure as
%       those used by mcx and mmc (which are both forward solvers).
%
%       The rbrunforward() function takes in most of the simulation
%       settings via the cfg struct input. This struct contains the following fields
%
%         node*: node list of the mesh
%         elem*: element list of the mesh
%         seg: labels/segmentation in the mesh (similar to mmc's cfg.elemprop)
%         prop*: wavelength-dependent optical properties
%         param: wavelength-independent physiological properties
%               (cfg.param.{hbo,hbr,water,lipids,aa3,scatamp,scatpow})
%         srctype*: source type [pencil]
%         srcpos*: source array locations
%         srcdir*: source directions
%         widesrc*: arrays of pattern source data (3D array: Nx*Ny*Npat)
%         dettype: detector type [pencil]
%         detpos*: detector array locations
%         detdir*: detector directions
%         bulk: bulk optical property
%         omega: frequency-domain system modulation (angular) frequency in rad
%
%       properties marked with "*" share nearly identical definitions as in
%       MMC (if exist, also in MCX).
%
%       The following extra inputs are derived from the above inputs -
%       they can be automatically populated by running rbmeshprep, or
%       manually added
%
%         face: mesh exterior surface triangles
%         area: areas of surface triangles
%         evol: element volume
%         nvol: nodal volume
%         reff: effective reflective coeff (computed by rbgetreff)
%         deldotdel: $\nabla\phi(r)\cdot\nabla\phi(r)$ - the quadratic term
%           of the FEM equation, (computed by rbdeldotdel)
%         {rows,cols,idxcount}: Redbird's FEM equation uses a Compressed
%           Sparse Column (CSC) format to store the system matrix A, these
%           are related to the sparse matrix columns returned by rbfemnz
%         idxsum: the cumsum() of cfg.idxcount, i.e. the starting index of
%           the non-zero terms for each node in the serialized vector
%
%     one can also pass on the cfg data structure used by mcxlab or mmclab
%     to rbrunforward. When cfg.nphoton is set, rbrunforward routes to a
%     Monte Carlo solver:
%       - cfg.node + cfg.elem  -> mmclab (tetrahedral mesh-based MMC)
%       - cfg.vol              -> mcxlab (voxel-grid MCX)
%     For the mmclab path, requesting the 6th output (Jext) additionally
%     computes the adjoint Jacobian via mmc's mesh-mode adjoint kernel
%     (requires cfg.detdir; auto-sets cfg.outputtype='adjoint_mua_d',
%     cfg.srcid=-1, cfg.basisorder=1). Forward-only calls (nargout<=5)
%     run as outputtype='fluence' and skip the Jacobian.
%     MWT/Helmholtz physics (cfg.helmholtz or cfg.bulkprop) is not
%     supported via the MC path.
%
% output:
%     detval: the values at the detector locations
%     phi: the full volumetric forward solution computed at all wavelengths
%     Amat: the left-hand-side matrices (a containers.Map object) at specified wavelengths
%     rhs: the right-hand-side vectors for all sources (independent of wavelengths)
%     sflag: solver flag (FEM solver), or 0 for the MC path
%     Jext: optional 6th output. Empty for the FEM path. For the MC path,
%           a struct with fields .mua (J_mua, [Nsd x Nn]) and .dcoeff (J_D,
%           [Nsd x Nn]); both fields wrap a containers.Map(wavelength) for
%           multi-spectral cfg.prop. Consumed by rbrunrecon to bypass rbjac.
%     param/value pairs: (optional) additional parameters
%          'solverflag': a cell array to be used as the optional parameters
%               for rbfemsolve (starting from parameter 'method'), for
%               example  rbrunforward(...,'solverflag',{'pcg',1e-10,200})
%               calls rbfemsolve(A,rhs,'pcg',1e-10,200) to solve forward
%               solutions
%
%
% license:
%     GPL version 3, see LICENSE_GPLv3.txt files for details
%
% -- this function is part of Redbird-m toolbox
%

opt = varargin2struct(varargin{:});

Jext = [];   % populated only by the MC path when nargout > 5

if (isfield(cfg, 'nphoton'))
    % MWT/Helmholtz is FEM-only physics
    if ((isfield(cfg, 'helmholtz') && cfg.helmholtz) || isfield(cfg, 'bulkprop'))
        error(['rbrunforward: MWT/Helmholtz forward (cfg.helmholtz / cfg.bulkprop) ' ...
               'cannot be combined with the Monte Carlo path. Remove cfg.nphoton.']);
    end

    % default time grid: only fill what the user left unset
    if (~isfield(cfg, 'tstart') || isempty(cfg.tstart))
        cfg.tstart = 0;
    end
    if (~isfield(cfg, 'tend') || isempty(cfg.tend))
        cfg.tend = 5e-9;
    end
    if (~isfield(cfg, 'tstep') || isempty(cfg.tstep))
        cfg.tstep = cfg.tend;     % single-gate (CW) by default
    end

    needJacobian = (nargout > 5);
    avgsize = jsonopt('avgsize', 1, opt);
    srcnum  = size(cfg.srcpos, 1);
    detnum  = size(cfg.detpos, 1);
    if (size(cfg.detpos, 2) == 3)
        cfg.detpos(:, 4) = avgsize;
    end

    if (isfield(cfg, 'node') && isfield(cfg, 'elem'))
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% mmclab path: tetrahedral mesh-based Monte Carlo
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        cfg.method     = 'elem';
        cfg.basisorder = 1;       % nodal fluence: required for mesh adjoint
        if (isfield(cfg, 'seg') && ~isfield(cfg, 'elemprop'))
            cfg.elemprop = cfg.seg;
        end

        if (needJacobian)
            % adjoint Jacobian path -> mmc returns J_mua and J_D
            if (~isfield(cfg, 'detdir') || isempty(cfg.detdir))
                if (jsonopt('autodetdir', 1, opt))
                    cfg.detdir = rbgetdetdir(cfg);
                    fprintf(1, '[rbrunforward] cfg.detdir auto-filled via rbgetdetdir (%d detectors)\n', size(cfg.detdir, 1));
                else
                    error(['rbrunforward: MC adjoint Jacobian requires cfg.detdir ' ...
                           '(Nd-by-4 inward detector normal + focal length).']);
                end
            end
            if (size(cfg.detdir, 2) < 4)
                cfg.detdir(:, 4) = 0;
            end
            if (isfield(cfg, 'omega') && cfg.omega > 0)
                error(['rbrunforward: mmc mesh-mode adjoint Jacobian does not yet ' ...
                       'support RF (cfg.omega > 0). Remove cfg.nphoton to use FEM.']);
            end
            cfg.srcid      = -1;
            cfg.outputtype = 'adjoint_mua_d';   % dual output: J_mua + J_D
        elseif (isfield(cfg, 'detdir') && ~isempty(cfg.detdir))
            % forward only, but user provided detdir: simulate detectors as
            % adjoint sources too (so phi has Ns+Nd slots) without computing J
            if (size(cfg.detdir, 2) < 4)
                cfg.detdir(:, 4) = 0;
            end
            cfg.srcid      = -2;
            cfg.outputtype = 'fluence';
        else
            cfg.outputtype = 'fluence';
        end

        % wavelength iteration (mmclab is single-wavelength per call)
        ismulti = isa(cfg.prop, 'containers.Map');
        if (ismulti)
            wavelengths = cfg.prop.keys;
            propAll     = cfg.prop;
            phi         = containers.Map();
            detval      = containers.Map();
            if (needJacobian)
                JmuaMap = containers.Map();
                JdMap   = containers.Map();
            end
        else
            wavelengths = {''};
        end

        [optodeloc, optodebary] = tsearchn(cfg.node(:, 1:3), cfg.elem(:, 1:4), cfg.detpos(:, 1:3));

        % Save the original multi-source cfg.srcpos / cfg.srcdir so the
        % per-source mmclab loop below can rotate one source at a time.
        % (mmclab takes cfg.srcpos as a single float4; multi-source redbird
        % configurations must be looped on the host side.)
        srcposAll = cfg.srcpos;
        srcdirAll = cfg.srcdir;

        for waveid = wavelengths
            wv = waveid{1};
            if (ismulti)
                cfg.prop = propAll(wv);
            end

            % Per-node optical-property mode for DOT reconstruction.
            % When cfg.prop has one row per forward-mesh node, route the
            % mua (and, for RF, musp) columns into cfg.nodemua / cfg.nodemusp
            % so mmc's per-node global-memory path is used (see mmc PR#3).
            % cfg.prop itself is reduced to a 2-row table
            % [outside; bulk-fallback] so the constant-memory gmed[] still
            % carries n and g for outside-mesh transitions; cfg.elemprop is
            % set to the single bulk tissue index everywhere.
            if (size(cfg.prop, 1) == size(cfg.node, 1))
                % keep as double - mmclab.cpp's parser reads via mxGetPr
                cfg.nodemua    = double(cfg.prop(:, 1));
                cfg.isnodalmua = 1;
                if (isfield(cfg, 'omega') && cfg.omega > 0)
                    cfg.nodemusp    = double(cfg.prop(:, 2));
                    cfg.isnodalmusp = 1;
                end
                bulk = mean(cfg.prop, 1);
                if (bulk(2) < 1e-3)
                    bulk(2) = 1e-3;
                end   % avoid mus=0 for path sampling
                cfg.prop = [0 0 1 1; bulk(1) bulk(2) 0 bulk(4)];
                cfg.elemprop = ones(size(cfg.elem, 1), 1);
            end

            phiwv    = zeros(size(cfg.node, 1), srcnum + detnum);
            detphiwv = nan(detnum, srcnum);
            if (needJacobian)
                Jmua_wv = zeros(srcnum * detnum, size(cfg.node, 1));
                Jd_wv   = zeros(srcnum * detnum, size(cfg.node, 1));
            end

            for s = 1:srcnum
                cfg.srcpos = srcposAll(s, :);
                if (size(srcdirAll, 1) > 1)
                    cfg.srcdir = srcdirAll(s, :);
                else
                    cfg.srcdir = srcdirAll;
                end
                cfg.srcdata = [];   % force mmclab to (re)build srcdata for this single source

                out = mmclab(cfg);
                if (needJacobian && ~isfield(out, 'jmua'))
                    fprintf(2, '[diag] s=%d: out fields = %s; outputtype=%s; size(detpos)=[%d,%d]; size(detdir)=[%d,%d]\n', ...
                            s, strjoin(fieldnames(out), ','), cfg.outputtype, ...
                            size(cfg.detpos, 1), size(cfg.detpos, 2), ...
                            size(cfg.detdir, 1), size(cfg.detdir, 2));
                end

                % flux.data shape this call: [nn, maxgate, 1+Nd] for multi-slot;
                % squeeze to drop maxgate=1 -> [nn, 1+Nd].
                phisingle = squeeze(out.data);
                if (~isvector(phisingle) && size(phisingle, 1) ~= size(cfg.node, 1) && size(phisingle, 2) == size(cfg.node, 1))
                    phisingle = phisingle.';
                end

                % forward-source fluence (slot 0) and detector-as-adjoint-source fluences (slots 1..Nd)
                phiwv(:, s) = phisingle(:, 1);
                phiwv(:, srcnum + 1:end) = phiwv(:, srcnum + 1:end) + phisingle(:, 2:end) / srcnum;

                % detphi at each detector position from this source
                for d = 1:detnum
                    if (~isnan(optodeloc(d)))
                        nodes_d = cfg.elem(optodeloc(d), 1:4);
                        detphiwv(d, s) = optodebary(d, :) * phisingle(nodes_d, 1);
                    end
                end

                if (needJacobian)
                    % out.jmua is [nn, 1*Nd] this call; map into rows (s-1)*Nd+1..s*Nd
                    jmua_s = double(out.jmua).';     % [Nd, nn]
                    jd_s   = double(out.jd).';
                    rowidx = (s - 1) * detnum + 1:s * detnum;
                    Jmua_wv(rowidx, :) = jmua_s;
                    Jd_wv(rowidx, :)   = jd_s;
                end
            end

            % restore multi-source cfg fields for any downstream caller
            cfg.srcpos = srcposAll;
            cfg.srcdir = srcdirAll;

            if (ismulti)
                phi(wv)    = phiwv;
                detval(wv) = detphiwv;
                if (needJacobian)
                    JmuaMap(wv) = Jmua_wv;
                    JdMap(wv)   = Jd_wv;
                end
            else
                phi    = phiwv;
                detval = detphiwv;
                if (needJacobian)
                    Jext = struct('mua', Jmua_wv, 'dcoeff', Jd_wv);
                end
            end
        end

        if (ismulti)
            cfg.prop = propAll;
            if (needJacobian)
                Jext = struct('mua', JmuaMap, 'dcoeff', JdMap);
            end
        end

        if (nargout > 2)
            Amat  = [];
        end
        if (nargout > 3)
            rhs   = [];
        end
        if (nargout > 4)
            sflag = 0;
        end

    elseif (isfield(cfg, 'vol'))
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% mcxlab path: voxel-grid Monte Carlo
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        cfg.srcid = -1;
        outdata   = cell(1, max(nargout, 1));
        [outdata{:}] = mcxlab(cfg);
        phi = squeeze(outdata{1}.data);
        detval = zeros(1, srcnum * detnum);
        for i = 1:size(phi, 4)
            detval(((i - 1) * detnum + 1):i * detnum) = rbvoxelmean(phi(:, :, :, i), cfg.detpos(:, 1:3), avgsize);
        end
        if (nargout > 2)
            Amat = outdata{3};
            if (nargout > 3)
                rhs = outdata{4};
                if (nargout > 4)
                    sflag = outdata{5};
                end
            end
        end
    else
        error('rbrunforward: cfg.nphoton requires either cfg.node+cfg.elem (mmclab) or cfg.vol (mcxlab)');
    end
    return
end

% time-domain DOT (Crank-Nicolson) is engaged when cfg.tstart, cfg.tstep,
% and cfg.tend are all defined; rbruntd handles the time loop.
if (isfield(cfg, 'tstart') && ~isempty(cfg.tstart) && ...
    isfield(cfg, 'tstep') && ~isempty(cfg.tstep) && ...
    isfield(cfg, 'tend') && ~isempty(cfg.tend))
    [detval, phi, Amat, rhs] = rbruntd(cfg, varargin{:});
    sflag = 0;
    return
end

rfcw = jsonopt('rfcw', 1, opt);

if (~isfield(cfg, 'deldotdel'))
    cfg.deldotdel = rbdeldotdel(cfg);
end

wavelengths = {'_'};   % non-empty, non-numeric placeholder for the single-wavelength case (kept internal; unwrapped before return)
if (isa(cfg.prop, 'containers.Map'))
    wavelengths = cfg.prop.keys;
end
if isfield(opt, 'sd') && ~isempty(opt.sd)
    sd = opt.sd;
else
    sd = rbsdmap(cfg);
end
if (~isa(sd, 'containers.Map'))
    sd = containers.Map(wavelengths, {sd});
end

Amat = containers.Map();
opt = varargin2struct(varargin{:});
solverflag = jsonopt('solverflag', {}, opt);

for md = rfcw
    detval(md).detphi = containers.Map();
    phi(md).phi = containers.Map();
end

sdtmp = cell2mat(sd.values');
srcnum = length(unique(sdtmp(:, 1)));

for waveid = wavelengths
    wv = waveid{1};

    sdwv = sd(wv);
    for md = rfcw

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%   Build RHS
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        rhs = rbfemrhs(cfg, sd, wv, md);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%   Build LHS
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        Amat(wv) = rbfemlhs(cfg, cfg.deldotdel, wv, md); % use native matlab code, 1 sec for 50k nodes

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%   Solve for solutions at all nodes: Amat*res=rhs
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % solverflag={'pcg',1e-12,200}; % if iterative pcg method is used
        [phi(md).phi(wv), sflag] = rbfemsolve(Amat(wv), rhs, solverflag{:});

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%   Extract detector readings from the solutions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        tempdetval = rbfemgetdet(phi(md).phi(wv), cfg, rhs); % or detval=rbfemgetdet(phi(wv), cfg, rhs);
        if size(sdwv, 2) < 4
            sdwv(:, 4) = ones(size(sdwv, 1), 1) .* md;
        end
        sdmd = sdwv(sdwv(:, 4) == md | sdwv(:, 4) == 3, :);
        detval(md).detphi(wv) = tempdetval(unique(sdmd(:, 2)) - srcnum, unique(sdmd(:, 1)));
    end
end

% if only a single wavelength is required, return regular arrays instead of a map
if (length(wavelengths) == 1)
    Amat = Amat(wavelengths{1});
    for md = rfcw
        phi(md).phi = phi(md).phi(wavelengths{1});
        detval(md).detphi = detval(md).detphi(wavelengths{1});
    end
end

if (length(rfcw) == 1)
    phi = phi(rfcw).phi;
    detval = detval(rfcw).detphi;
end
