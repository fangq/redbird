function Jscat = rbjacscat(Jd, dcoeff, scatpow, wv, lref, suffix)
%
% Jscat=rbjacscat(Jd, dcoeff, scatpow, wv, lref, suffix)
%
% Building the Jacobian matrices for scattering-amplitude and
% scattering-power from Jacobian of the diffusion coeff
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     Jd: the Jacobian for diffusion coefficient, as a containers.Map
%     dcoeff: the current diffusion coeff at each node
%     scatpow: the current scattering power at each node
%     wv: wavelength list, if not given, Jd must be a containers.Map
%     lref: reference wavelength (in nm) used in the inverse power law
%           musp = scatamp * (wavelen/lref)^(-scatpow). Default: 1e9 (legacy
%           formula). Pass 500 for the 500 nm-normalized convention.
%     suffix: suffix appended to the output field names, default ''.
%             Pass '500' so the fields are named scatamp500/scatpow500
%             when using the 500 nm-normalized convention.
%
% output:
%     Jscat: the Jacobian in a struct as Jscat.{scatamp,scatpow}
%            (or Jscat.{scatamp500,scatpow500} when suffix='500')
%
% license:
%     GPL version 3, see LICENSE_GPLv3.txt files for details
%
% -- this function is part of Redbird-m toolbox
%

if (nargin < 4)
    wv = keys(Jd);
end
if (nargin < 5 || isempty(lref))
    lref = 1e9;
end
if (nargin < 6 || isempty(suffix))
    suffix = '';
end

ampname = ['scatamp' suffix];
powname = ['scatpow' suffix];

% Jscat=[J(scatamp),J(scatpow)]
if isa(Jd, 'containers.Map')
    Jscat = struct(ampname, [], powname, []);
    for i = 1:length(wv)
        Jscat.(ampname) = [Jscat.(ampname); rbjacscatamp(Jd(wv{i}), dcoeff(wv{i}), str2double(wv{i}), scatpow, lref)];
        Jscat.(powname) = [Jscat.(powname); rbjacscatpow(Jd(wv{i}), dcoeff(wv{i}), str2double(wv{i}), lref)];
    end
elseif isstruct(Jd)
    Jscat = struct(ampname, {[], []}, powname, {[], []});
    for ii = 1:2
        for ll = 1:length(wv)
            Jscat(ii).(ampname) = [Jscat(ii).(ampname); rbjacscatamp(Jd(ii).J(wv{ll}), dcoeff(wv{ll}), str2double(wv{ll}), scatpow, lref)];
            Jscat(ii).(powname) = [Jscat(ii).(powname); rbjacscatpow(Jd(ii).J(wv{ll}), dcoeff(wv{ll}), str2double(wv{ll}), lref)];
        end
    end
end
% Jscat.scatamp=zeros(size(cell2mat(Jd.values'),1),size(cell2mat(Jd.values'),2));
% Jscat.scatpow=Jscat.scatamp;
