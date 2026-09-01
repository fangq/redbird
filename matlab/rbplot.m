function varargout = rbplot(cfg, data, varargin)
%
% h=rbplot(cfg, data, varargin)
%
% Plot the forward solution or Jacobian defined by the 2nd parameter data
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     cfg: the simulation settings stored as a redbird data structure
%     data: the data to be plotted - the row number can be either #node or #elem
%
% output:
%     h: the handles of the plotted graphics objects
%
% license:
%     GPL version 3, see LICENSE_GPLv3.txt files for details
%
% -- this function is part of Redbird toolbox
%

if (nargin == 1)
    cfg.detpos(:, 4) = 1;
    mcxpreview(cfg);
    return
end
if (size(data, 2) == size(cfg.node, 1) || size(data, 2) == size(cfg.elem, 1))
    [varargout{1:nargout}] = rbplotjacobian(cfg, data, varargin{:});
elseif (size(data, 1) == size(cfg.node, 1) || size(data, 1) == size(cfg.elem, 1))
    [varargout{1:nargout}] = rbplotforward(cfg, data, varargin{:});
else
    error('data row or column counts must match the node or elem number');
end
