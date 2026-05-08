function [sd, dist] = rbsdmap(cfg, varargin)
%
% [sd,dist]=rbsdmap(cfg,maxdist)
% [sd,dist]=rbsdmap(cfg,'maxdist',10,'excludesrc',[1,2],'excludedet',[1],'wavesrc',...,'wavedet',...)
%
% Create a source-detector mapping table (sd) by setting a maximum separation distance
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     cfg: the redbird data structure
%     maxdist: the maximum source detector separation distance (in mm)
%
% output:
%     sd: the source detector mapping table
%     dist: the full table of all source detector separations (#row=#sources, #col=#detectors)
%
% license:
%     GPL version 3, see LICENSE_GPLv3.txt files for details
%
% -- this function is part of Redbird-m toolbox
%

if ((~isfield(cfg, 'srcpos') && ~isfield(cfg, 'widesrc')) || (isempty(cfg.srcpos) && ~isfield(cfg, 'widesrc')) || (~isfield(cfg, 'detpos') && ~isfield(cfg, 'widedet')) || (isempty(cfg.detpos) && ~isfield(cfg, 'widedet')))
    error('you must define at least 1 source and 1 detector in order to use this function');
end

if (length(varargin) == 1)
    maxdist = varargin{1};
elseif (isempty(varargin))
    maxdist = inf;
elseif (length(varargin) > 1)
    opt = varargin2struct(varargin{:});
    maxdist = jsonopt('maxdist', inf, opt);
end

if ~exist('opt', 'var')
    opt = struct();
end

srcnum = size(cfg.srcpos, 1);
detnum = size(cfg.detpos, 1);
if isfield(cfg, 'widesrc')
    widesrcnum = size(cfg.widesrc, 1);
    if isfield(cfg, 'excludewfsrc')
        badwfsrc = cfg.excludewfsrc;
    elseif isfield(opt, 'excludewfsrc')
        badwfsrc = jsonopt('excludewfsrc', [], opt);
    else
        badwfsrc = [];
    end
    badwfsrc = unique([badwfsrc find((sum(cfg.widesrc, 2) == 0) | isnan(sum(cfg.widesrc, 2)))']);
else
    widesrcnum = 0;
    badwfsrc = [];
end
if isfield(cfg, 'widedet')
    widedetnum = size(cfg.widedet, 1);
    if isfield(cfg, 'excludewfdet')
        badwfdet = cfg.excludewfdet;
    elseif isfield(opt, 'excludewfdet')
        badwfdet = jsonopt('excludewfdet', [], opt);
    else
        badwfdet = [];
    end
    badwfdet = unique([badwfdet find((sum(cfg.widedet, 2) == 0) | isnan(sum(cfg.widedet, 2)))']);
else
    widedetnum = 0;
    badwfdet = [];
end

if isfield(cfg, 'excludesrc')
    badsrc = cfg.excludesrc;
else
    badsrc = jsonopt('excludesrc', [], opt);
end
if isfield(cfg, 'excludedet')
    baddet = cfg.excludedet;
else
    baddet = jsonopt('excludedet', [], opt);
end

if (isfield(cfg, 'srcpos') && ((size(cfg.srcpos, 2) == size(cfg.face, 1)) || (size(cfg.srcpos, 2) == size(cfg.node, 1))))
    dist = zeros(srcnum, detnum);
else
    if isfield(cfg, 'widesrc')
        widesrc = cfg.widesrc;
    else
        widesrc = [];
    end
    if isfield(cfg, 'widedet')
        widedet = cfg.widedet;
    else
        widedet = [];
    end
    dist = rbgetdistance(cfg.srcpos, cfg.detpos, badsrc, baddet, widesrc, widedet, cfg)';
end

src = 1:srcnum;
widesrc = [1:widesrcnum] + srcnum;
% normalize to 1xN row vectors: setdiff returns a 0x1 empty when the input
% is empty in some Octave versions (>=10), which breaks horizontal
% concatenation with non-empty row vectors below.
goodsrc = reshape(sort(setdiff(1:srcnum, badsrc)), 1, []);
goodwfsrc = reshape(sort(setdiff(1:widesrcnum, badwfsrc)), 1, []);
if ~isempty(goodwfsrc)
    goodwfsrc = goodwfsrc + srcnum;
end
det = 1:detnum;
widedet = [1:widedetnum] + detnum;
gooddet = reshape(sort(setdiff(1:detnum, baddet)), 1, []);
goodwfdet = reshape(sort(setdiff(1:widedetnum, badwfdet)), 1, []);
if ~isempty(goodwfdet)
    goodwfdet = goodwfdet + detnum;
end

if (isfield(cfg, 'prop') && isa(cfg.prop, 'containers.Map'))
    wavelengths = cfg.prop.keys;
    sd = containers.Map();
    [ss, dd] = meshgrid([src widesrc], srcnum + widesrcnum + [det widedet]);

    for wv = wavelengths
        wid = wv{1};
        if (isfield(cfg, 'wfsrcmapping'))
            wfsrcmap = cfg.wfsrcmapping(wid);
        end
        if (isfield(cfg, 'wfdetmapping'))
            wfdetmap = cfg.wfdetmapping(wid);
        end

        if ((isfield(cfg, 'wavesrc') && ~isempty(cfg.wavesrc)) || ...
            (isfield(cfg, 'wavedet') && ~isempty(cfg.wavedet)))
            if (~isempty(cfg.wavesrc(wid)))
                wavesrc = cfg.wavesrc(wid);
                if exist('wfsrcmap', 'var')
                    wavesrc = rbremapsrc(wavesrc, wfsrcmap, srcnum);
                end
                goodwavesrc = intersect([goodsrc goodwfsrc], wavesrc);
            else
                wavesrc = goodsrc;
            end
            if (~isempty(cfg.wavedet(wid)))
                wavedet  = cfg.wavedet(wid);
                if exist('wfdetmap', 'var')
                    wavedet = rbremapsrc(wavedet, wfdetmap, detnum);
                end
                goodwavedet = intersect([gooddet goodwfdet], wavedet);
            else
                wavedet = gooddet;
            end
            [ss, dd] = meshgrid(wavesrc, srcnum + widesrcnum + wavedet);
        else
            goodwavesrc = [goodsrc goodwfsrc];
            goodwavedet = [gooddet goodwfdet];
        end
        sdwv = [ss(:), dd(:)];
        %         if(nargin<2 || (size(cfg.srcpos,2) == size(cfg.face,1)))
        if (~isinf(maxdist))
            [s2, d2] = meshgrid(1:srcnum, (1:detnum) + srcnum + widesrcnum);
            sd2 = [s2(:) d2(:)];
            [tf, idx] = ismember(sdwv, sd2, 'rows');
            idx = idx(find(idx));
            sdwv(:, 3) = 1;
            sdwv(:, 3) = (reshape(dist(unique(sdwv(:, 2)) - (srcnum + widesrcnum), unique(sdwv(:, 1))), [], 1) < maxdist);
            % sdwv(idx,3)=(dist(:)<maxdist);
        else
            sdwv(:, 3) = 0;
            sdwv(ismember(sdwv(:, 1), goodwavesrc) & ismember(sdwv(:, 2), goodwavedet + srcnum + widesrcnum), 3) = 1;
        end

        if (isfield(cfg, 'rfcw'))
            modes = cfg.rfcw.src.keys;
            sdwv(:, 4) = zeros(size(sdwv, 1), 1);
            for md = modes
                mid = md{1};
                mdSRC = cfg.rfcw.src(mid);
                if (exist('wfsrcmap', 'var') && any(ismember(wfsrcmap(:, 1), mdSRC)))
                    mdSRC = rbremapsrc(mdSRC, wfsrcmap, srcnum);
                end
                mdDET = cfg.rfcw.det(mid);
                if (exist('wfdetmap', 'var') && any(ismember(wfdetmap(:, 1), mdDET)))
                    mdDET = rbremapsrc(mdDET, wfdetmap, detnum);
                end
                mdDET = mdDET + srcnum + widesrcnum;

                mdChan = (ismember(sdwv(:, 1), mdSRC) & ismember(sdwv(:, 2), mdDET));
                if strcmp(mid, 'RF')
                    sdwv(mdChan, 4) = sdwv(mdChan, 4) + 1;
                elseif strcmp(mid, 'CW')
                    sdwv(mdChan, 4) = sdwv(mdChan, 4) + 2;
                end
            end
            sdwv(sdwv(:, 4) == 0, :) = [];
            %             sdwv((sdwv(:,3) == 0 | sdwv(:,4) == 0),:) = [];
        else
            sdwv((sdwv(:, 3) == 0), :) = [];
        end
        sd(wid) = sdwv;
    end
else
    [ss, dd] = meshgrid([goodsrc goodwfsrc + srcnum], [gooddet goodwfdet + detnum] + (srcnum + widesrcnum));
    sd = [ss(:), dd(:)];
    if (nargin < 2 || (size(cfg.srcpos, 2) == size(cfg.face, 1)))
        sd(:, 3) = 1;
    else
        sd(:, 3) = (dist(:) < maxdist);
    end
end

function remap = rbremapsrc(sourcelist, wflist, srcnum)
if (isempty(intersect(sourcelist, wflist)) && ~isempty(wflist))
    remap = sourcelist;
else
    [idA, idB] = ismember(wflist(:, 1), sourcelist);
    remap = sourcelist;
    remap(idB(find(idB))) = [];
    wf = wflist(idA, :);
    for ii = 1:size(wf, 1)
        remap = [remap srcnum + (wf(ii, 2):wf(ii, 3))];
    end
end
