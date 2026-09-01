%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Redbird - A Diffusion Solver for Diffuse Optical Tomography,
%      Copyright Qianqina Fang, 2018
%
% This example shows explicitly the detailed steps of running a forward
% simulation. One can call rbrun or rbrunforward as one-liner alternatives
%
% This file is part of Redbird URL:http://mcx.sf.net/mmc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (~exist('rbrun', 'file'))
    addpath(fullfile(pwd, '../matlab'));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   prepare simulation input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear cfg;

[cfg.node, cfg.face, cfg.elem] = meshabox([40 0 0], [160, 120, 60], 10);

% [cfg.node, cfg.face, cfg.elem]=meshabox([0 0 0],[60 60 30],3);
nn = size(cfg.node, 1);
cfg.seg = ones(size(cfg.elem, 1), 1);

[xi, yi] = meshgrid(60:20:140, 20:20:100);
cfg.srcpos = [xi(:), yi(:), zeros(numel(yi), 1)];
cfg.detpos = [xi(:), yi(:), 60 * ones(numel(yi), 1)];
cfg.srcdir = [0 0 1];
cfg.detdir = [0 0 -1];

cfg.prop = [
            0 0 1 1
            0.008 1 0 1.37
            0.016 1 0 1.37
           ];

cfg.omega = 2 * pi * 70e6;
cfg.omega = 0;

tic;
cfg = rbmeshprep(cfg);
fprintf('preparing mesh ... \t%f seconds\n', toc);

% save config.mat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Build LHS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic;
deldotdel = rbdeldotdel(cfg);
Amat = rbfemlhs(cfg, deldotdel); % use native matlab code, 1 sec for 50k nodes
fprintf('build LHS using native code ... \t%f seconds\n', toc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Build RHS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[rhs, loc, bary] = rbfemrhs(cfg);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Solve for solutions at all freenodes: Afree*sol=rhs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic;
fprintf(1, 'solving for the solution ...\n');
% [phi,sflag]=rbfemsolve(Amat,rhs,'pcg',1e-8,200);
phi = rbfemsolve(Amat, rhs);
fprintf('solving forward solutions ... \t%f seconds\n', toc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Extract detector readings from the solutions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

detval = rbfemgetdet(phi, cfg, loc, bary); % or detval=rbfemgetdet(phi, cfg, rhs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Analytical solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sid = 13;

srcloc = cfg.srcpos(sid, 1:3);
detloc = cfg.node;

phicw = cwdiffusion(cfg.prop(2, 1), cfg.prop(2, 2) * (1 - cfg.prop(2, 3)), cfg.reff, srcloc, detloc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Visualization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
subplot(221);
plotmesh([cfg.node, log10(abs(phi(1:size(cfg.node, 1), sid)))], cfg.elem, 'y=30', 'facecolor', 'interp', 'linestyle', 'none');
cl = get(gca, 'clim');
set(gca, 'xlim', [60, 140]);
set(gca, 'zlim', [0 60]);
view([0 1 0]);
colorbar;

subplot(222);
plotmesh([cfg.node, log10(abs(phicw(1:size(cfg.node, 1), 1)))], cfg.elem, 'y=30', 'facecolor', 'interp', 'linestyle', 'none');
view([0 1 0]);
set(gca, 'xlim', [60, 140]);
set(gca, 'zlim', [0 60]);
set(gca, 'clim', cl);
colorbar;

dd = log10(abs(phi(1:size(cfg.node, 1), sid))) - log10(abs(phicw(1:size(cfg.node, 1), 1)));
subplot(223);
plotmesh([cfg.node, dd], cfg.elem, 'y=30', 'facecolor', 'interp', 'linestyle', 'none');
view([0 1 0]);
set(gca, 'xlim', [60, 140]);
set(gca, 'zlim', [0 60]);
colorbar;

subplot(224);
plotmesh([cfg.node, dd], cfg.elem, 'y=30', 'facecolor', 'interp', 'linestyle', 'none');
hist(dd(:), 100);

%% test add-noise function

dist = rbgetdistance(cfg.srcpos, cfg.detpos);
plot(dist(:), log10(abs(detval(:))), '.');
% snrshot is referenced to the strongest channel (max|detval|~2e-7 here);
% this 43 dB matches the 110 dB used before rbaddnoise was renormalized
newdata = rbaddnoise(detval, 43, 40);
hold on;
plot(dist(:), log10(abs(newdata(:))), 'r.');
hold off;
