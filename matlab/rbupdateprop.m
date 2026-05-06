function prop = rbupdateprop(cfg, wv)
%
% prop=rbupdateprop(cfg)
%
% Update the direct material properties (cfg.prop - optical or EM
% properties used by the forward solver) using multispectral properties
% (cfg.param - physiological parameters) for all wavelengths
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     cfg: the redbird data structure
%
% output:
%     prop: the updated material property data to replace cfg.prop
%
% license:
%     BSD or GPL version 3, see LICENSE_GPLv3.txt files for details
%
% -- this function is part of Redbird-m toolbox
%

% for single-wavelength
if (~isfield(cfg, 'param') && (~isa(cfg.prop, 'containers.Map') || length(cfg.prop) == 1))
    prop = cfg.prop;
    return
end

% fieldnames(cfg.param) provides chromorphore species, 2nd input wv or
% keys(cfg.prop) provides the wavelength list; cfg.prop is a containers.Map
% object and cfg.param is a struct.

if (~isfield(cfg, 'prop') || ~isa(cfg.prop, 'containers.Map'))
    error('input cfg must be a struct and must have subfield names "prop" and "param"');
end

if (nargin < 2)
    wv = keys(cfg.prop);
end

prop = containers.Map();
params = cfg.param;

for i = 1:length(wv)
    wavelen = wv{i};
    if ~isfield(params, 'water')
        params.water = 0.23;
    end
    if ~isfield(params, 'lipids')
        params.lipids = 0.58;
    end
    types = intersect(fieldnames(params), {'hbo', 'hbr', 'water', 'lipids', 'aa3'});
    if (isempty(types))
        error('specified parameters are not supported');
    end
    extin = rbextinction(str2double(wavelen), types);
    mua = zeros(size(params.(types{1})));
    for j = 1:length(types)
        mua = mua + extin(j) * params.(types{j});
    end
    if (isfield(cfg.param, 'scatamp500') && isfield(cfg.param, 'scatpow500'))
        musp = (cfg.param.('scatamp500') .* ((str2double(wavelen) ./ 500).^(-cfg.param.('scatpow500'))));
    elseif (isfield(cfg.param, 'scatamp') && isfield(cfg.param, 'scatpow'))
        musp = (cfg.param.('scatamp') .* ((str2double(wavelen) .* 1e-9).^(-cfg.param.('scatpow'))));
    end
    segprop = cfg.prop(wavelen);
    if (length(mua) < min([size(cfg.node, 1), size(cfg.elem, 1)])) % label-based properties
        segprop(length(mua) + 2:end, :) = [];
        segprop(2:end, 1) = mua(:);
        if (exist('musp', 'var'))
            segprop(2:end, 2) = musp(:);
            segprop(2:end, 3) = 0;
        end
        prop(wavelen) = segprop;
    else % mua and musp are defined on the node or elements
        if (exist('musp', 'var'))
            prop(wavelen) = [mua(:) musp(:) zeros(size(musp(:))) segprop(2, 4) * ones(size(musp(:)))];
        else
            segprop = repmat(segprop(2, 2:end), length(mua), 1);
            prop(wavelen) = [mua(:) segprop];
        end
    end
end
