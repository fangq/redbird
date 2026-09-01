function nodeval = rbelem2node(elem, elemval, nodelen)
%
% nodeval=rbelem2node(elem,elemval,nodelen)
%
% Interpolaing solutions from element-based values to node-based values
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     elem: the mesh element list
%     elemval: values defined per element, a matrix of size Ne x Nv, where
%         Ne - number of elements
%         Nv - number of values per element
%     nodelen: total number of nodes
%
% output:
%     nodeval: the interpolated values at each node
%
% license:
%     GPL version 3, see LICENSE_GPLv3.txt files for details
%
% -- this function is part of Redbird toolbox
%

if (isstruct(elem))
    nodelen = size(elem.node, 1);
    elem = elem.elem;
end

nelem = size(elem);
nval = size(elemval);

if (nval(2) == nelem(1))
    nodeval = zeros(nval(1), nodelen);
    for j = 1:nelem(2)
        nodeval(:, elem(:, j)) = nodeval(:, elem(:, j)) + elemval(:, elem(:, j));
    end
else
    nodeval = zeros(nodelen, nval(1));
    for j = 1:nelem(2)
        nodeval(elem(:, j), :) = nodeval(elem(:, j), :) + elemval(elem(:, j), :);
    end
end
nodeval = nodeval * 0.25;
