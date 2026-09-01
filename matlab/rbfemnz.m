function [rows, cols, connnum] = rbfemnz(elem, nn)
%
% [rows,cols,connnum]=rbfemnz(elem,nn)
%
% Return the indices of all non-zero elements in a sparse matrix (for FEM)
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     elem:  element table of a mesh
%     nn:    number of all nodes
%
% output:
%     rows: output, the row index of the non-zero off-diagonal elements
%     cols: output, the column index of the non-zero off-diagonal elements
%     connnum: output, total non-zero off-diagonal element number
%
%     newdata: noise contaminated data
%
% license:
%     GPL version 3, see LICENSE_GPLv3.txt files for details
%
% -- this function is part of Redbird toolbox
%

[conn, connnum, count] = meshconn(elem, nn);

rows = zeros(1, count);
cols = zeros(1, count);
count = 0;
for i = 1:nn
    rows(count + 1:count + connnum(i)) = i;
    cols(count + 1:count + connnum(i)) = conn{i};
    count = count + connnum(i);
end
