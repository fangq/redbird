function varargout = rbplotforward(cfg, phi, idx, varargin)
%
% h=rbplotforward(cfg, phi, idx, 'param1',value1, 'param2',value2,...)
%
% Plot one column of the forward solution phi
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     cfg: the simulation settings stored as a redbird data structure
%     phi: the forward solution matrix
%     idx: the index of the column to be plotted
%     options: additional plotting options in the form of 'parame', value pairs.
%
% output:
%     h: the handles of the plotted graphics objects
%
% license:
%     GPL version 3, see LICENSE_GPLv3.txt files for details
%
% -- this function is part of Redbird toolbox
%

if (nargin < 3)
    idx = 1;
end

len = length(idx);

hd = zeros(len, 1);

for i = 1:len
    subplot(1, len, i);
    hd(i) = plotmesh([cfg.node, phi(:, idx(i))], cfg.elem, varargin{:});
end

if (nargout > 0)
    varargout{1} = hd;
end
