function [Jmua_n, Jmua_e, Jd_n, Jd_e] = rbjac(sd, phi, deldotdel, felem, evol, varargin)
%
% [Jmua_n, Jmua_e, Jd_n, Jd_e]=rbjac(sd, phi, deldotdel, felem, evol, isnodal)
%
% Building the Jacobian matrices using adjoint method and native matlab code
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     sd: the source-detector mapping table
%     phi: the forward solutions at all sources
%     deldotdel: grad*phi dot product with grad phi, computed as part of the computation
%     felem: forward mesh element list
%     evol: forward mesh element volume
%     iselem (optional): 1 if elem-based Jacobian; 0 (default) if is nodal-based
%
% output:
%     Jmua_n: the nodal Jacobian for absorption coeff. mua
%     Jmua_e: the element-wise Jacobian for absorption coeff. mua
%     Jd_n: (optional) the nodal Jacobian for diffusion coeff D
%     Jd_e: (optional) the element-wise Jacobian for diffusion coeff D
%
% license:
%     GPL version 3, see LICENSE_GPLv3.txt files for details
%
% -- this function is part of Redbird-m toolbox
%

if (nargin < 5 || isempty(sd) || isempty(phi) || isempty(deldotdel) || isempty(evol))
    error('you must give at least 5 inputs and they must not be empty');
end
if (isempty(varargin))
    iselem = 0;
    rfcw = 1;
elseif (length(varargin) == 1)
    iselem = varargin{1};
    rfcw = 1;
else
    opt = varargin2struct(varargin{:});
    rfcw = jsonopt('rfcw', 1, opt);
    iselem = jsonopt('iselem', 0, opt);
end
nelem = size(felem, 1);

wavelengths = {'_'};

if isstruct(phi)
    if (isa(phi(1).phi, 'containers.Map'))
        wavelengths = phi(1).phi.keys;
    else
        for ii = 1:length(phi)
            phi(ii).phi = containers.Map({'_'}, {phi(ii).phi});
        end
    end
else
    if (isa(phi, 'containers.Map'))
        wavelengths = phi.keys;
        phitemp = phi;
        clear phi;
        phi(rfcw).phi = phitemp;
        clear phitemp;
    else
        phitemp = phi;
        clear phi;
        phi(rfcw).phi = containers.Map({'_'}, {phitemp});
        clear phitemp;
    end
end
mode = 1:length(phi);

for ii = mode
    Jmua_n(ii).J = containers.Map();
    Jmua_e(ii).J = containers.Map();
    if (nargout > 2)
        Jd_n(ii).J = containers.Map();
        Jd_e(ii).J = containers.Map();
    end
end

idx = [1 1 1 2 2 3 2 3 4 3 4 4
       2 3 4 3 4 4 1 1 1 2 2 3
       2 3 4 6 7 9 2 3 4 6 7 9];

for waveid = wavelengths
    wv = waveid{1};
    if (isa(sd, 'containers.Map'))
        sdwv = sd(wv);
    else
        sdwv = sd;
    end
    for md = rfcw
        if (size(sdwv, 2) < 4)
            sdwv(:, 4) = md;
        end
        sdmd = sdwv(((sdwv(:, 3) == 1) & (sdwv(:, 4) == md | sdwv(:, 4) == 3)), :);
        phiwv = phi(md).phi(wv);
        Jmua_node = zeros(size(sdmd, 1), size(phiwv, 1));
        Jmua_elem = zeros(size(sdmd, 1), nelem);
        if (nargout > 2)
            Jd_node = zeros(size(sdmd, 1), size(phiwv, 1));
            Jd_elem = zeros(size(sdmd, 1), nelem);
        end
        for i = 1:nelem
            phidotphi1 = phiwv(felem(i, 1:4), sdmd(:, 1)) .* phiwv(felem(i, 1:4), sdmd(:, 2));
            phidotphi2 = phiwv(felem(i, idx(1, :)), sdmd(:, 1)) .* phiwv(felem(i, idx(2, :)), sdmd(:, 2));
            Jmua_elem(:, i) = -(sum(phidotphi1, 1) + sum(phidotphi2, 1) * 0.5) * (0.1 * evol(i));
            if (nargout > 2)
                Jcol = phidotphi1.' * deldotdel(i, [1 5 8 10]).';
                Jcol = Jcol + phidotphi2.' * deldotdel(i, idx(3, :)).';
                Jd_elem(:, i) = -Jcol;
            end
            for j = 1:4
                Jmua_node(:, felem(i, j)) = Jmua_node(:, felem(i, j)) + Jmua_elem(:, i);
            end
            if (nargout > 2)
                for j = 1:4
                    Jd_node(:, felem(i, j)) = Jd_node(:, felem(i, j)) - Jcol;
                end
            end
        end
        Jmua_node = Jmua_node * 0.25;
        Jmua_n(md).J(wv) = Jmua_node;
        Jmua_e(md).J(wv) = Jmua_elem;
        if (nargout > 2)
            Jd_node = Jd_node * 0.25;
            Jd_n(md).J(wv) = Jd_node;
            Jd_e(md).J(wv) = Jd_elem;
        end
    end
end

% if only a single wavelength is required, return regular arrays instead of a map
if (length(wavelengths) == 1)
    for md = mode
        Jmua_n(md).J = Jmua_n(md).J(wavelengths{1});
        Jmua_e(md).J = Jmua_e(md).J(wavelengths{1});
        if (nargout > 2)
            Jd_n(md).J = Jd_n(md).J(wavelengths{1});
            Jd_e(md).J = Jd_e(md).J(wavelengths{1});
        end
    end
end

if (length(rfcw) == 1)      %% (length(mode) == 1)
    Jmua_n = Jmua_n(rfcw).J;    %% Jmua_n(1).J;
    Jmua_e = Jmua_e(rfcw).J;    %% Jmua_e(1).J;
    if (nargout > 2)
        Jd_n = Jd_n(rfcw).J;    %% Jd_n(1).J;
        Jd_e = Jd_e(rfcw).J;    %% Jd_e(1).J;
    end
end

if (iselem && nargout < 3)
    Jmua_n = Jmua_e;
    if nargout > 2
        Jmua_e = Jd_e;
    end
end
