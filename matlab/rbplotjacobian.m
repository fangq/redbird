function varargout = rbplotjacobian(cfg, jac, idx, varargin)
%
% h=rbplotjacobian(cfg, jac, idx, 'param1',value1, 'param2',value2,...)
%
% Plot one row of the Jacobian matrix
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     cfg: the simulation settings stored as a redbird data structure
%     jac: the Jacobian matrix
%     idx: the index of the row to be plotted
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
    if (size(jac, 2) == size(cfg.node, 1))
        hd(i) = plotmesh([cfg.node, jac(idx(i), :)'], cfg.elem, varargin{:});
    elseif (size(jac, 2) == size(cfg.elem, 1))
        hd(i) = plotmesh(cfg.node, [cfg.elem(:, 1:4), jac(idx(i), :)'], varargin{:});
    end
end

if (nargout > 0)
    varargout{1} = hd;
end
