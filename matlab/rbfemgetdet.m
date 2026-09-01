function [detval, goodidx] = rbfemgetdet(phi, cfg, optodeloc, optodebary)
%
% [detval, goodidx]=rbfemgetdet(phi, cfg, rhs)
%    or
% [detval, goodidx]=rbfemgetdet(phi, cfg, optodeloc, optodebary)
%
% Retrieving measurement data at detectors as a (N_det by N_src) matrix
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     phi: the forward solution obtained by rbforward or rbfemsolve
%     cfg: the redbird simulation data structure
%     rhs: the RHS matrix returned by rbfemrhs (only use this when the
%          source or detector contains widefield sources)
%     optodeloc: the optode enclosing element ID returned by rbgetoptodes
%     optodebary: the optode barycentric coordinates returned by rbgetoptodes
%
% output:
%     detval: #det x #src array, denoting the measurement data
%     goodidx: the list of "good optodes" - for point optodes, the non-NaN
%        optodeloc indices; for widefield sources- all src/det are included
%
% license:
%     GPL version 3, see LICENSE_GPLv3.txt files for details
%
% -- this function is part of Redbird toolbox
%

if (((~isfield(cfg, 'srcpos') || isempty(cfg.srcpos)) && ~isfield(cfg, 'widesrc')) || ...
    ((~isfield(cfg, 'detpos') || isempty(cfg.detpos)) && ~isfield(cfg, 'widedet')))
    detval = [];
    goodidx = [];
    return
end

srcnum = size(cfg.srcpos, 1);
detnum = size(cfg.detpos, 1);
if isfield(cfg, 'widesrc')
    wfsrcnum = size(cfg.widesrc, 1);
else
    wfsrcnum = 0;
end
if isfield(cfg, 'widedet')
    wfdetnum = size(cfg.widedet, 1);
else
    wfdetnum = 0;
end

goodsrc = find(~isnan(optodeloc(1:srcnum + wfsrcnum)));
goodidx = find(~isnan(optodeloc(srcnum + wfsrcnum + 1:srcnum + wfsrcnum + detnum + wfdetnum)));
detval = zeros(length(goodidx), length(goodsrc));
gooddetval = zeros(srcnum, detnum);

if (nargin == 3)
    %     detval=optodeloc(:,srcnum+1:srcnum+detnum)'*phi(:,1:srcnum);
    [tf, goodsrc] = find(sum(optodeloc(:, 1:srcnum + wfsrcnum)));
    goodsrc = unique(goodsrc);
    [tf, goodidx] = find(sum(optodeloc(:, srcnum + wfsrcnum + 1:srcnum + wfsrcnum + detnum + wfdetnum)));
    goodidx = unique(goodidx);
    % use non-conjugating transpose: the widefield detector RHS is complex for
    % MWT (scaled by -j*omega*mu0); a conjugate transpose here would flip the
    % sign of the measurement relative to the adjoint-based Jacobian in rbjac
    detval = optodeloc(:, [1:detnum + wfdetnum] + srcnum + wfsrcnum).' * phi(:, [1:srcnum + wfsrcnum]);
    %     detval=optodeloc(:,goodidx+srcnum+wfsrcnum)'*phi(:,goodsrc);
elseif (isempty(goodidx) && size(cfg.detpos, 2) == size(cfg.node, 1)) % wide-field det
    for i = 1:srcnum
        for j = 1:detnum
            detval(j, i) = sum(phi(:, i) .* cfg.detpos(j, :)');
        end
    end
else
    if (~isempty(goodidx))
        for i = goodidx'
            %         if(~isnan(optodeloc(goodidx(i))))
            gooddetval(i, goodsrc) = sum(phi(cfg.elem(optodeloc(srcnum + i), :), goodsrc) .* repmat(optodebary(srcnum + i, :)', 1, length(goodsrc)), 1);
            %         end
        end
    end
    detval = gooddetval(goodidx, goodsrc);
end
