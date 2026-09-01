function varargout = rbfemmatrix(varargin)
% [Adiag, Aoff, deldotdel]=rbfemmatrix(cfg)
%   or
% [Jmua, Jd]=rbfemmatrix(cfg, sd, phi, deldotdel, isnodal)
% [Jmua, Jd]=rbfemmatrix(cfg, sd, phi, deldotdel, isnodal, mapid, mapweight, rnodenum, relem)
%
% Mex-based core function to build FEM forward matrix and Jacobian matrix
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     sd: the source-detector mapping table
%     phi: the forward solutions at all sources
%     deldotdel: grad*phi dot product with grad phi, computed as part of
%           the computation
%     isnodal: 0 for element-based Jacobian and 1 for node-based Jacobian
%     mapid: the list of element indices of the reconstruction mesh where
%           each forward mesh node is enclosed
%     mapweight: the barycentric coordinates of the forward mesh nodes
%           inside the reconstruction mesh elements
%     rnodenum: total node number of the reconstruction mesh
%     relem: reconstruction mesh element list
%
% output:
%     Adiag: diagonal portion of the FEM forward matrix
%     Aoff: non-zero off-diagonal portion of the FEM forward matrix
%     Jmua: the nodal or elementary Jacobian for absorption coeff. mua
%     Jd: (optional) the node or elementary Jacobian for diffusion coeff D
%
% license:
%     GPL version 3, see LICENSE_GPLv3.txt files for details
%
% -- this function is part of Redbird toolbox
%
