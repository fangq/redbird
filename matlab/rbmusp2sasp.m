function [sa, sp] = rbmusp2sasp(musp, lambda)
%
% [sa,sp]=rbmusp2sasp(musp,lambda)
%
% Computing scattering amplitude and scattering power based on mus' at two
% wavelengths, they are connected by lambda=sa*(lambda/500 nm)^(-sp)
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     musp: reduced scattering coeff in 1/mm, a vector of length 2
%     lambda: wavelengths in nm, a vector of length 2
%
% output:
%     sa: scattering amplitude
%     sp: scattering power
%
% license:
%     GPL version 3, see LICENSE_GPLv3.txt files for details
%
% -- this function is part of Redbird toolbox
%

lambda = lambda / 500;
sp = log(musp(1) / musp(2)) / log(lambda(2) / lambda(1));
sa = 0.5 * (musp(1) / lambda(1)^(-sp) + musp(2) / lambda(2)^(-sp));
