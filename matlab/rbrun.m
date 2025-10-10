function varargout = rbrun(cfg, recon, detphi0, varargin)
%
% detphi=rbrun(cfg)
% newrecon=rbrun(cfg,recon,detphi0)
% [newrecon,resid,newcfg,...]=rbrun(cfg,recon,detphi0,sd)
%   or
% [newrecon,resid,newcfg,...]=rbrun(cfg,recon,detphi0,sd,'param1',value1,'param2',value2,...)
%
% Master script to run streamlined forward, inversion with various
% reconstruction modes.
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     cfg: the forward data structure defining the forward mesh in a
%          reconstruction, if cfg is the only input, rbrun is the same as
%          rbrunforward
%     recon: recon is a struct, it can have the following subfields
%               recon.node: the node list of the reconstruction mesh
%               recon.elem: the element list of the reconstruction mesh
%               recon.mapid: a vector of length the forward mesh node,
%                     specifying the recon mesh element ID that encloses
%                     each of the forward mesh node; NaN if the forward
%                     mesh is outside of the recon mesh
%               recon.mapweight: a Nn-by-4 array, Nn is the forward mesh
%                     node number. Each row defines the 4 barycentric
%                     coordinates for each of the forward mesh node; all
%                     NaN if the forward mesh is outside of the recon mesh
%               recon.lambda: Tikhonov regularization parameter, if not
%                     present, use 0.05
%               recon.bulk: a struct defining the fall-back values of the
%                     bulk optical properties, can have the following
%                     subfields:
%                       hbo/hbr/scatamp/scatpow: bulk HbO/HbR
%                           concentrations and scattering amp/power, if
%                           defined, it will be used to
%                           initialize/overwrite recon.param then cfg.param
%                           then cfg.prop
%                       mua/musp/dcoeff: bulk absorption coeff (1/mm),
%                           reduced scattering coeff (1/mm) and diffusion
%                           coeff. (mm), if defined, it will be used to
%                           initialize overwrite recon.prop then cfg.prop
%                       n: refractive index, if not defined, use 1.37
%
%                       recon.bulk must be defined if some forward nodes
%                       are outside of the recon mesh
%               recon.param: initial guess of multi-spectral parameters
%                     (wavelength independent) such as chromophore
%                     concentrations and scattering amp/power;
%                     this must be defined, but can be an empty struct;
%                     rbrun needs the subfield names from this struct to
%                     determine which parameters to use
%               recon.prop: initial guess of wavelength-dependent
%                     properties, always have 4 columns -
%                           [mua (1/mm), musp (1/mm), g, n]
%                     anisotropy g is often set to 0; n is set to 1.37 if
%                     recon.bulk.n is not defined.
%
%                     recon.prop must be a matrix (for single-wavelength
%                     recon) or a containers.Map object (for multi-spectral
%                     recon). it must be defined if multi-spectral recon is
%                     used to provide the desired wavelengths as
%                     recon.prop.keys
%     detphi0: measurement data to be fitted, can be either a matrix or a
%              containers.Map (multi-wavelength); if detphi0 is a struct
%              with subfields like node/elem, it is treated as a forward
%              simulation data structure like cfg, and is used to run the
%              forward solution to obtain the simulated "measurement"
%
% output:
%     newrecon,resid,newcfg...: if only cfg is used as input, the output
%              parameters are the same as in rbrunforward; if recon is
%              given, the output sequence is the same as that of rbrunrecon
%
% license:
%     GPL version 3, see LICENSE_GPLv3.txt files for details
%
% -- this function is part of Redbird-m toolbox
%

mode = 'image';
opt = struct;
if (nargin == 4 || (~isempty(varargin) && (~ischar(varargin{1}) || bitand(length(varargin), 1) == 1 && ~ischar(varargin{1}))))
    sd = varargin{1};
    if (length(varargin) == 2 && ~ischar(varargin{1}) && ischar(varargin{2}))
        mode = varargin{2};
    end
    if (length(varargin) >= 3)
        opt = varargin2struct(varargin{2:end});
    end
else
    if (length(varargin) == 2 && ~ischar(varargin{1}) && ischar(varargin{2}))
        sd = varargin{1};
        mode = varargin{2};
    else
        opt = varargin2struct(varargin{:});
    end
end

mode = jsonopt('mode', mode, opt);
method = jsonopt('method', 'normal', opt);
prior = jsonopt('prior', '', opt);
if isfield(opt, 'sd')
    sd = opt.sd;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Run forward for the heterogeneous domain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (nargin == 1)
    [varargout{1:nargout}] = rbrunforward(cfg);
    return
elseif (nargin == 2)
    if (~isfield(recon, 'detphi0'))
        error('you must give detphi0 as input');
    else
        detphi0 = recon.detphi0;
    end
end
if (isstruct(detphi0) && isfield(detphi0, 'srcpos'))
    detphi0 = rbrunforward(detphi0);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Reset the domain to a homogeneous medium for recon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch mode
    case {'bulk', 'seg', 'image'}
        if (strcmp(mode, 'bulk'))
            if (isfield(recon, 'node'))
                recon.seg = ones(size(recon.node, 1), 1);
            else
                cfg.seg = ones(size(cfg.node, 1), 1);
            end
            maxseg = 1;
        elseif (strcmp(mode, 'image'))
            if (~isempty(prior) && isfield(recon, 'seg') && ~isfield(opt, 'lmat'))
                opt.lmat = rbprior(recon.seg, prior, opt);
            end
            if (isfield(recon, 'seg'))
                recon = rmfield(recon, 'seg');
            end
            if (isfield(cfg, 'seg'))
                cfg = rmfield(cfg, 'seg');
            end
            if (isfield(recon, 'node'))
                maxseg = size(recon.node, 1);
            else
                maxseg = size(cfg.node, 1);
            end
        else
            if (isfield(recon, 'node'))
                maxseg = max(cfg.seg);
            else
                maxseg = max(cfg.seg);
            end
        end
        if (isfield(recon, 'param') && isstruct(recon.param) && isfield(recon, 'bulk'))
            types = fieldnames(recon.param);
            for ch = 1:length(types)
                recon.param.(types{ch}) = recon.bulk.(types{ch}) * ones(maxseg, 1);
            end
            if (strcmp(mode, 'image'))
                for ch = 1:length(types)
                    cfg.param.(types{ch}) = recon.bulk.(types{ch}) * ones(size(cfg.node, 1), 1);
                end
            end
            % cfg.prop will be updated inside rbrunrecon;
        end
        if ((~isfield(recon, 'param') && isfield(recon, 'prop') && isempty(recon.prop)) || ...
            (~isfield(recon, 'param') && ~isfield(recon, 'prop')))
            nref = 1.37;
            if (isfield(recon.bulk, 'n'))
                nref = recon.bulk.n;
            end
            if (isfield(recon.bulk, 'musp'))
                recon.prop = repmat([recon.bulk.mua, recon.bulk.musp, 0, nref], maxseg, 1);
            else
                recon.prop = repmat([recon.bulk.mua, 1 / (3 * recon.bulk.dcoeff), 0, nref], maxseg, 1);
            end
            if (strcmp(mode, 'image'))
                if (isfield(recon.bulk, 'musp'))
                    cfg.prop = repmat([recon.bulk.mua, recon.bulk.musp, 0, nref], size(cfg.node, 1), 1);
                else
                    cfg.prop = repmat([recon.bulk.mua, 1 / (3 * recon.bulk.dcoeff), 0, nref], size(cfg.node, 1), 1);
                end
            end
            if (strcmp(mode, 'image') == 0)
                recon.prop = [0 0 1 1; recon.prop];
            end
        end
end

if (isfield(recon, 'seg') && ~isempty(recon.seg))
    if (numel(recon.seg) == 1 || length(unique(recon.seg)) == 1)
        recon.seg = ones(size(recon.node, 1), 1);
        fprintf('[rbrun]: run bulk estimation ...\n');
    elseif (isvector(recon.seg))
        fprintf('[rbrun]: run segmented tissue property estimation ...\n');
    elseif (min(size(recon.seg)) >= 2)
        if (strcmp(mode, 'fuzzy'))
            fprintf('[rbrun]: run fuzzy-segmentation based property estimation ...\n');
        elseif (strcmp(mode, 'prior'))
            fprintf('[rbrun]: run compositional-prior guided reconstruction ...\n');
        elseif (strcmp(mode, 'roi'))
            fprintf('[rbrun]: run roi-based compositional-prior guided reconstruction ...\n');
        end
    end
end

if (~exist('sd', 'var'))
    sd = rbsdmap(cfg);
end

[varargout{1:nargout}] = rbrunrecon(cfg, recon, detphi0, sd, 'maxiter', 10, opt);
