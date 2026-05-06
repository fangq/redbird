function [Jscatpow, dDdscatpow] = rbjacscatpow(Jd, dcoeff, wavelen, lref)
%
% [Jscatpow,dDdscatpow]=rbjacscatpow(Jd, dcoeff, wavelen, lref)
%
% Create the Jacobian matrix for the scattering power using the Jacobian
% of the diffusion coeff (D)
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     Jd: the jacobian of the diffusion coeff.
%     dcoeff: the diffusion coefficient values at each node
%     wavelen: the list of wavelengths (in nm)
%     lref: reference wavelength (in nm) used in the inverse power law
%           musp = scatamp * (wavelen/lref)^(-scatpow). Default: 1e9 (legacy
%           formula). Pass 500 for the 500 nm-normalized convention.
%
% output:
%     Jscatpow: the Jacobian of the scattering power parameter
%     dDdscatpow: partial derivative of D - diffusion coeff - to the scat power
%
% license:
%     GPL version 3, see LICENSE_GPLv3.txt files for details
%
% -- this function is part of Redbird-m toolbox
%

if (nargin < 4 || isempty(lref))
    lref = 1e9;
end

dDdscatpow = dcoeff .* log(wavelen ./ lref);

Jscatpow = Jd .* dDdscatpow;
