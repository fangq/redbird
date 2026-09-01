function [deldotdel, delphi] = rbdeldotdel(cfg)
%
% [deldotdel, delphi]=rbdeldotdel(cfg)
%
% Compute deldotdel=<grad(phi_i).grad(phi_j)>, where <> means spatial integration
% inside elements, "." means dot-product, grad() means gradience, phi means
% linear basis function in a tetrahedron. For a linear function phi,
% grad(phi) is a constant across the element.
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     cfg: the redbird data structure
%
% output:
%     deldotdel: precomputed deldotdel=<grad(phi_i).grad(phi_j)>. For each
%         element, deldotdel is a 4x4 symmetric matrix - to store this data
%         efficiently; we only store the upper triangule of the matrix per
%         element - this gives a Ne x 10 matrix (Ne is the number of tets);
%         for each row: [d11, d12, d13, d14, d22, d23, d24, d33, d34, d44]
%     delphi: gradient of the basis functions (grad(phi_i)) in each
%         tetrahedral element; a 3x4 matrix for each element with a total
%         dimension is 3 x 4 x Ne, where
%            3 - gradient direction, for x, y and z
%            4 - which basis functions - for node 1, 2, 3 and 4
%            Ne - number of element
%
% license:
%     GPL version 3, see LICENSE_GPLv3.txt files for details
%
% -- this function is part of Redbird toolbox
%

no = cfg.node;
el = cfg.elem(:, 1:4);

no = reshape(no(el', :)', 3, 4, size(el, 1));

delphi = zeros(3, 4, size(el, 1));

col = [4 2 3 2
       3 1 4 3
       2 4 1 4
       1 3 2 1];

for coord = 1:3
    idx = 1:3;
    idx(coord) = [];
    for i = 1:4
        delphi(coord, i, :) = squeeze(((no(idx(1), col(i, 1), :) - no(idx(1), col(i, 2), :)) .* (no(idx(2), col(i, 3), :) - no(idx(2), col(i, 4), :)) - ...
                                       (no(idx(1), col(i, 3), :) - no(idx(1), col(i, 4), :)) .* (no(idx(2), col(i, 1), :) - no(idx(2), col(i, 2), :)))) ./ (cfg.evol(:) * 6);
    end
end

deldotdel = zeros(size(el, 1), 10);
count = 1;

for i = 1:4
    for j = i:4
        deldotdel(:, count) = sum(squeeze(delphi(:, i, :) .* delphi(:, j, :)), 1);
        count = count + 1;
    end
end
deldotdel = deldotdel .* repmat(cfg.evol(:), 1, 10);
