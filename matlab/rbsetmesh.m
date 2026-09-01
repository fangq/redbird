function cfg = rbsetmesh(cfg0, node, elem, prop, propidx)
%
% cfg=rbsetmesh(cfg0,node,elem,prop,propidx)
%
% Associate a new tetrahedral mesh to a redbird simulation structure
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     cfg0: the initial redbird data structure
%     node: the node coordinate list of the new mesh, an Nn-by-3 floating point matrix
%     elem: the elem index list of the new mesh, an Ne-by-4 integer matrix
%     prop: the direct property list
%     propidx: the property mapping list
%
% output:
%     cfg: the updated redbird data structure
%
% license:
%     GPL version 3, see LICENSE_GPLv3.txt files for details
%
% -- this function is part of Redbird toolbox
%

names = fieldnames(cfg0);
names = intersect(names, {'face', 'area', 'evol', 'deldotdel', 'isreoriented', 'nvol', 'cols', 'idxsum', 'musp0'});

cfg0.node = node;
cfg0.elem = elem;
if (nargin > 3)
    cfg0.prop = prop;
end
if (nargin > 4)
    cfg0.seg = propidx;
end

if (~isempty(names))
    cfg = rmfield(cfg0, names);
else
    cfg = cfg0;
end
cfg = rbmeshprep(cfg);
