function cfgkey = rbgetcfg(cfgs, key)
%
% cfgkey=rbgetcfg(cfgs,wavelength)
%
% Creating wavelength-specific simulation structure from a containers.Map
% based multi-spectral simulation structure
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     cfgs: multi-spectral simulation structure, if cfgs does not contain
%         multi-spectral information, return value is the same as input
%     wavelength: a string or a number indicating the wavelength to be
%         extracted; this wavelength must exist in the cfgs struct
%
% output:
%     cfgkey: the single-wavelength cfg structure at the given wavelength
%
% license:
%     GPL version 3, see LICENSE_GPLv3.txt files for details
%
% -- this function is part of Redbird toolbox
%

if (~ischar(key))
    key = num2str(key);
end

cfgkey = structfun(@(x) getkey(x, key), cfgs);

function newx = getkey(x, key)
newx = x;
if (isa(x, 'containers.Map'))
    newx = x(key);
end
