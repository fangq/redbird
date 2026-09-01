function [Anew, allkeys] = rbmatflat(Amat, weight)
%
% [Anew, allkeys]=rbmatflat(Amat)
%
% Flatten a containers.Map object (in the case of multi-wavelength cases)
% into a rectangular 2D matrix
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     Amat: the LHS of the matrix equation
%     weight: the weight mutiplying to each sub block before concatenation
%
% output:
%     Anew: the flattened LHS matrix
%     allkeys: the wavelength labels in the Amat input if it is a containers.Map, or []
%
% license:
%     GPL version 3, see LICENSE_GPLv3.txt files for details
%
% -- this function is part of Redbird toolbox
%

allkeys = [];
if (~isstruct(Amat) && ~isa(Amat, 'containers.Map'))
    Anew = Amat;
    return
end

Anew = [];
if (isa(Amat, 'containers.Map')) % maps are vertically concatenated (wavelength)
    allkeys = Amat.keys;
    if (nargin < 2)
        weight = ones(length(allkeys), 1);
    end
    for i = 1:length(allkeys)
        Anew = [Anew; Amat(allkeys{i}) * weight(i)];
    end
elseif (isstruct(Amat)) % structs are horizontally concatenated (chromorphores)
    allkeys = fieldnames(Amat);
    rfcw = length(Amat);
    Anew = struct('J', cell(size(Amat)));
    if (nargin < 2)
        weight = ones(length(allkeys), 1);
    end
    for i = 1:length(allkeys)
        for j = 1:rfcw
            Anew(j).J = [Anew(j).J, Amat(j).(allkeys{i}) * weight(i)];
        end
    end
    if (length(Anew) == 1)
        Anew = Anew(1).J;
    end
elseif (iscell(Amat)) % cells are combined via cell2mat
    Anew = cell2mat(Amat);
end
