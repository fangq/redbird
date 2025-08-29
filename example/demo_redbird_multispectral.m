%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Redbird - A Diffusion Solver for Diffuse Optical Tomography,
%      Copyright Qianqina Fang, 2018
%
% In this example, we show the most basic usage of Redbird.
%
% This file is part of Redbird URL:http://mcx.sf.net/mmc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (~exist('rbrun', 'file'))
    addpath(fullfile(pwd, '../matlab'));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   prepare simulation input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear cfg cfg0 recon;

s0 = [70, 50, 20]; % center of the inclusion (in mm)
rs = 5;            % radius of the sphere (in mm)

[nobbx, fcbbx] = meshabox([40 0 0], [160, 120, 60], 10);
[nosp, fcsp] = meshasphere(s0, rs, 1);
[no, fc] = mergemesh(nobbx, fcbbx, nosp, fcsp);

[cfg0.node, cfg0.elem] = s2m(no, fc(:, 1:3), 1, 40, 'tetgen', [41 1 1; s0]);

nn = size(cfg0.node, 1);
cfg0.seg = cfg0.elem(:, 5);
cfg0.srcdir = [0 0 1];

[xi, yi] = meshgrid(60:20:140, 20:20:100);
cfg0.srcpos = [xi(:), yi(:), zeros(numel(yi), 1)];
cfg0.detpos = [xi(:), yi(:), 60 * ones(numel(yi), 1)];
cfg0.detdir = [0 0 -1];

% cfg.param defines the wavelength-independent optical/physiological properties
% cfg.param accepts: hbo (in uM), hbr (in uM), water (in 0-1 volume fraction), lipids (in 0-1 volume fraction) - these are used to compute mua at any wavelength
% cfg.param also accepts: scatamp, scatpow - these are used to compute mus' at any wavelength via mus'=scatamp*lambda(in m)^scatpow

cfg0.param = struct;
cfg0.param.hbo = [15 30];
cfg0.param.hbr = [4  8];

% cfg.prop is the wavelength-specific optical properties, it is a volatile variable to be updated when running each wavelength
% if both prop and param are defined, param will be used to update/ovewrite prop; the function to map param to prop is rbupdateprop()
% one should still define an initial version of prop so that it can be updated

cfg0.prop = containers.Map();

% the below two lines serve the following purposes
% 1. the keys in cfg.prop provides the wavelength list, as strings of numbers (can have digits)
% 2. it has N+1 rows, for tissue type 0, and type 1-N, just like MCX/MMC; the first row is for type 0
% 3. it defines refractive index (n) for each tissue type
% 4. anisotropy g must be 0 in redbird because it solves the DE
% here, because param.hbo and param.hbr are defined, the mua of the below cfg.prop will be updated (the initial values do not matter here)
% if param.scatamp and param.scatpow are also defined, the mus will also be updated (therefore, the initial values will be overwritten)
cfg0.prop('690') = [0 0 1 1; 0   1 0 1.37; 0 1 0 1.37];
cfg0.prop('830') = [0 0 1 1; 0 0.8 0 1.37; 0 0.8 0 1.37];

wavelengths = cfg0.prop.keys;

cfg0.omega = 2 * pi * 70e6;
cfg0.omega = 0;

cfg = cfg0;

cfg0 = rbmeshprep(cfg0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   run forward for all wavelengths
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

detphi0 = rbrun(cfg0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   run reconstruction using the forward data, setup dual-mesh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[node, face, elem] = meshabox([40 0 0], [160, 120, 60], 10);
clear face;

cfg = rbsetmesh(cfg, node, elem, cfg.prop, ones(size(node, 1), 1));

[recon.node, face, recon.elem] = meshabox([40 0 0], [160, 120, 60], 25);
clear face;
[recon.mapid, recon.mapweight] = tsearchn(recon.node, recon.elem, cfg.node);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   run bulk fitting first
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sd = rbsdmap(cfg);
recon.bulk = struct('hbo', 8, 'hbr', 2); % Required: this gives initial guesses
recon.param = struct('hbo', 8, 'hbr', 2); % Required: this defines chromophores
recon.prop = containers.Map({'690', '830'}, {[], []}); % Required: for wavelengths
[newrecon, resid] = rbrun(cfg, recon, detphi0, sd, 'mode', 'bulk', 'lambda', 1e-3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   take the fitted bulk and set it for full image recon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

recon.bulk.hbo = newrecon.param.hbo;
recon.bulk.hbr = newrecon.param.hbr;
[newrecon, resid, newcfg] = rbrun(cfg, recon, detphi0, sd, 'mode', 'image', 'lambda', 1e-3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Plotting results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plotmesh([newrecon.node, newrecon.param.hbo(:)], newrecon.elem, 'z=20', 'facecolor', 'interp', 'linestyle', 'none');
hold on;
plotmesh([newrecon.node, newrecon.param.hbo(:)], newrecon.elem, 'x=70', 'facecolor', 'interp', 'linestyle', 'none');
view(3);
