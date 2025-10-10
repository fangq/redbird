function [detval, phi, Amat, rhs, sflag] = rbrunforward(cfg, varargin)
%
% [detval, phi]=rbrunforward(cfg)
%    or
% [detval, phi, Amat, rhs]=rbrunforward(cfg,'param1',value1,...)
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
%     to rbrunforward. If rbrunforward detects cfg.nphoton and
%     cfg.node/cfg.elem in the input, it calls mmclab and run forward
%     simulations on both the source (cfg.srcpos) and detectors
%     (cfg.detpos); if it sees cfg.nphoton and cfg.vol, then, rbrunforward
%     calls mcxlab to run the forward simulation. mcxlab/mmclab can not
%     handle frequency-domain simulations
%
% output:
%     detval: the values at the detector locations
%     phi: the full volumetric forward solution computed at all wavelengths
%     Amat: the left-hand-side matrices (a containers.Map object) at specified wavelengths
%     rhs: the right-hand-side vectors for all sources (independent of wavelengths)
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

if (isfield(cfg, 'nphoton'))
    if (~isfield(cfg, 'tend'))
        cfg.tend = 5e-9;
    end
    cfg.tstep = cfg.tend;
    cfg.outputtype = 'fluence';
    outdata = cell(1, nargout);
    srcnum = size(cfg.srcpos, 1);
    detnum = size(cfg.detpos, 1);
    avgsize = jsonopt('avgsize', 1, opt);
    detval = zeros(1, srcnum * detnum);
    srcdetnum = srcnum + detnum;
    if (nargout > 2)
        Amat = cell(1, srcdetnum);
    end
    if (size(cfg.detpos, 2) == 3)
        cfg.detpos(:, 4) = avgsize;
    end
    if (isfield(cfg, 'node') && isfield(cfg, 'elem'))
        cfg.method = 'elem';
        if (isfield(cfg, 'seg') && ~isfield(cfg, 'elemprop'))
            cfg.elemprop = cfg.seg;
        end
        srcpos = [cfg.srcpos(:, 1:3); cfg.detpos(:, 1:3)];
        [optodeloc, optodebary] = tsearchn(cfg.node(:, 1:3), cfg.elem(:, 1:4), cfg.detpos(:, 1:3));
        phi = zeros(size(cfg.node, 1), srcdetnum);
        for i = 1:srcdetnum
            cfg.srcpos = srcpos(i, :);
            [outdata{:}] = mmclab(cfg);
            phi(:, i) = outdata{1}.data;
            if (nargout > 2)
                Amat{i} = outdata{2};
            end
        end
        detval = sum(phi(cfg.elem(optodeloc, 1:4)) .* optodebary, 2);
    elseif (isfield(cfg, 'vol'))
        cfg.srcid = -1;
        [outdata{:}] = mcxlab(cfg);
        phi = squeeze(outdata{1}.data);
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
        error('input cfg is not supported by redbird, mcxlab and mmclab');
    end
    return
end

rfcw = jsonopt('rfcw', 1, opt);

if (~isfield(cfg, 'deldotdel'))
    cfg.deldotdel = rbdeldotdel(cfg);
end

wavelengths = {''};
if (isa(cfg.prop, 'containers.Map'))
    wavelengths = cfg.prop.keys;
end
sd = jsonopt('sd', containers.Map(wavelengths, cell(1, length(wavelengths))), opt);
if (isempty(sd(wavelengths{1})))
    sd = rbsdmap(cfg);
    if (~isa(sd, 'containers.Map'))
        sd = containers.Map(wavelengths, {sd});
    end
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
