function [Jmua_node, Jd_node] = rbjacnode(Jmua_elem, Jd_elem, elem, nodelen)
%
% [Jmua_node, Jd_node]=rbjacnode(Jmua_elem, Jd_elem, elem, nodelen)
%
% Converting from element-based Jacobians to node-based Jacobians
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     Jmua_elem: the element-wised Jacobian of the absorption coeff.
%     Jd_elem: the element-wised Jacobian of the diffusion coeff.
%     elem: the mesh element list
%     nodelen: the total number of the mesh nodes
%
% output:
%     Jmua_node: the nodal Jacobian of the absorption coeff.
%     Jd_node: the nodal Jacobian of the diffusion coeff.
%
% license:
%     GPL version 3, see LICENSE_GPLv3.txt files for details
%
% -- this function is part of Redbird toolbox
%

if (nargout > nargin)
    error('output number can not exceed input number');
end

Jmua_node = rbelem2node(elem, Jmua_elem, nodelen);

if (nargout > 1 && nargin > 1)
    Jd_node = rbelem2node(elem, Jd_elem, nodelen);
end
