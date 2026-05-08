function Jmua = rbjacmuafast(sd, phi, nvol, elem)
%
% Jmua=rbjacmuafast(sd, phi, nvol, elem)
%
% computing the nodal Jacobian matrix for mua using an approximated
% method as shownn in Qianqian's Thesis (nodal-adjoint formulation)
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     sd: the source-detector mapping table
%     phi: the forward solutions at all sources
%     nvol: the nodal volume computed by nodevolume from iso2mesh
%     elem: the mesh element list
%
% output:
%     Jmua: the approximated Jacobian for mua using nodal-adjoint formulation
%
% license:
%     GPL version 3, see LICENSE_GPLv3.txt files for details
%
% -- this function is part of Redbird-m toolbox
%

if (nargin < 3 || isempty(sd) || isempty(phi) || isempty(nvol))
    error('you must give at least the first 3 inputs and they must not be empty');
end

wavelengths = {'_'};

if (isa(phi, 'containers.Map'))
    wavelengths = phi.keys;
else
    phi = containers.Map({'_'}, {phi});
end

Jmua = containers.Map();

for waveid = wavelengths
    wv = waveid{1};
    phiwv = phi(wv);
    if (isa(sd, 'containers.Map'))
        sdwv = sd(wv);
    else
        sdwv = sd;
    end
    if (size(phiwv, 1) == length(nvol))
        Ja = zeros(size(sdwv, 1), size(phiwv, 1));
        for i = 1:size(sdwv, 1)
            Ja(i, :) = phiwv(:, sdwv(i, 1)) .* phiwv(:, sdwv(i, 2)) .* nvol(:);
        end
    elseif (nargin > 3 && size(nvol, 1) == size(elem, 1))
        Ja = zeros(size(elem, 1), size(sdwv, 1));
        for i = 1:size(sdwv, 1)
            for j = 1:4
                Ja(:, i) = Ja(:, i) + phiwv(elem(:, j), sdwv(i, 1)) .* phiwv(elem(:, j), sdwv(i, 2)) .* nvol(:);
            end
        end
        Ja = Ja.' * 0.25;
    else
        error('the row number of phi must be the same as the length of nvol or elem is needed');
    end
    Jmua(wv) = -Ja; % increasing mua, decreasing phi
end

% if only a single wavelength is required, return regular arrays instead of a map
if (length(wavelengths) == 1)
    Jmua = Jmua(wavelengths{1});
end
