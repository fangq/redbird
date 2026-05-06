function [Jscatamp, dDdscatamp] = rbjacscatamp(Jd, dcoeff, wavelen, scatpow, lref)
%
% [Jscatamp,dDdscatamp]=rbjacscatamp(Jd, dcoeff, wavelen, scatpow, lref)
%
% Create the Jacobian matrix for the scattering amplitude using the Jacobian
% of the diffusion coeff (D)
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     Jd: the jacobian of the diffusion coeff.
%     dcoeff: the diffusion coefficient values at each node
%     wavelen: the list of wavelengths (in nm)
%     scatpow: the scattering power of the current estimate of scattering power
%     lref: reference wavelength (in nm) used in the inverse power law
%           musp = scatamp * (wavelen/lref)^(-scatpow). Default: 1e9 (legacy
%           formula musp = scatamp*(wavelen-in-meters)^(-scatpow)). Pass 500
%           to use the 500 nm-normalized convention.
%
% output:
%     Jscatamp: the Jacobian of the scattering amplitude parameter
%     dDdscatamp: partial derivative of D - diffusion coeff - to the scat amplitude
%
% license:
%     GPL version 3, see LICENSE_GPLv3.txt files for details
%
% -- this function is part of Redbird-m toolbox
%

if (nargin < 5 || isempty(lref))
    lref = 1e9;
end

dDdscatamp = -3 .* dcoeff .* dcoeff .* (wavelen ./ lref).^(-scatpow.');

Jscatamp = Jd .* dDdscatamp;
