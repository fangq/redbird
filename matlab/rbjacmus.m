function Jmus = rbjacmus(Jd, mus, g)
%
% Jmus=rbjacmus(Jd, musp, g)
%
% Converting from the Jacobian of diffusion coeff to the Jacobian of mus'
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     Jd: the Jacobian of the diffusion coeff.
%     musp: the reduced scattering coeff
%     g: the anisotropy
%
% output:
%     Jmus: the nodal Jacobian of the absorption coeff.
%
% license:
%     GPL version 3, see LICENSE_GPLv3.txt files for details
%
% -- this function is part of Redbird toolbox
%

if (nargin < 3)
    g = 0;
end
Jmus = -Jd * (1 / (3 * mus * mus * (1 - g)));
