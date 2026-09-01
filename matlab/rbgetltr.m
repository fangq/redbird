function ltr = rbgetltr(cfg, wavelength)
%
% ltr=rbgetltr(cfg,wavelength)
%
% Compute the transport mean free path (l_tr=1/mu_tr) in mm in a medium
% where mu_tr=mua+musp is the transport coefficient, mua is the absorption
% coeff and musp=mus*(1-g) is the reduced scattering coeff, mus is the
% scattering coeff and g is the anisotropy
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     cfg: the forward simulation data structure
%     wavelength (optional): if cfg.prop is a containers.Map for
%          multispectral simulations, wavelength specifies which
%          wavelength, it can be a string or an integer.
%
% output:
%     ltr: transport mean free path in mm
%
% license:
%     GPL version 3, see LICENSE_GPLv3.txt files for details
%
% -- this function is part of Redbird toolbox
%

bkprop = rbgetbulk(cfg);
if (isa(bkprop, 'containers.Map'))
    if (nargin == 1)
        wavelength = bkprop.keys;
        wavelength = wavelength{1};
    end
    bkprop = bkprop(wavelength);
end
ltr = 1 / (bkprop(1) + bkprop(2) * (1 - bkprop(3)));
