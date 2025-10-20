clear cfg

disp('== mesh generation (iso2mesh) ...')
%# Create a box-like homogeneous domain
boxsize = [60, 50, 40]; %# domain size
trisize = 4; %# max triangule size on the surface
maxvol = 4; %# max tet element volume

tic
[nbox1, fbox1] = meshabox([0,0,0], boxsize, trisize);        %# create a big box
[nbox2, fbox2] = meshabox([10,10,10], [30,30,30], trisize);  %# create a box inclusion
[nbox1, fbox1] = removeisolatednode(nbox1, fbox1(:,1:3));    %# clean the surface
[nbox2, fbox2] = removeisolatednode(nbox2, fbox2(:,1:3));

[no1, fc1] = mergemesh(nbox1, fbox1, nbox2, fbox2);          %# combine the two non-intersecting surfaces
regionseed = [1, 1, 1        %# a seed point inside region 1 (large box)
              11,11,11]      %# a seed point inside region 2 (small box)
[cfg.node, cfg.elem] = s2m(no1,fc1,1,maxvol, 'tetgen', regionseed);  %# generate tetrahedral mesh - outer box: label 1, inner box: label 2

cfg.seg = cfg.elem(:,5);    %# cfg.seg is similar to cfg.elemprop in mmclab, but also supports node-based labels
cfg.elem(:,5)=[];
toc

plotmesh(cfg.node, [cfg.elem, cfg.seg], 'y>20')

disp('== define settings in cfg ...')

%# Creating forward simulation data structure cfg
%# properties:[mua(1/mm),mus(1/mm),g, n]  % if both mus/g are given, mus'=mus/(1-g) will be used for diffusion, usually set g to 0
cfg.prop  =  [0,        0,        1, 1           %# cfg.prop is the same as mcx/mmc, row-1 is for label 0, row-2 for label 1 etc
              0.006,    0.8,      0, 1.37        %# label 1 (background domain)
              0.02,     1,      0, 1.37];      %# label 2 (inclusion)

%# Default redbird source is a pencil beam, same as mcx/mmc
cfg.srcpos = [25, 25, 0];    %# redbird srcpos can have multiple rows
cfg.srcdir = [0, 0, 1];          %# srcdir determines how sources are sunken into the mesh

%# Redbird detector positions are point-like, directly sampling the output fluence; it is different from the disk-like shape in mcx/mmc
cfg.detpos = [35, 25, max(cfg.node(:,3))];  %# redbird detpos can have multiple rows
cfg.detdir = [0, 0, -1];          %# redbird automatically computes the adjoint solutions, treating detectors as sources

%# redbird cfg does not need cfg.{nphoton,tstart,tend,tstep,elemprop}, which are required in mmclab

disp('== preprocessing domain ...')
%# calling rbmeshprep() populates other needed mesh data, such as cfg.{evol,nvol,face,area,reff,deldotdel,cols,idxsum}
%# this is similar to `cfg = mmclab(cfg, 'prep')`

tic
cfg = rbmeshprep(cfg);
toc

disp('== run simulation ...')
%# Run forward simulation (you can also call rbrun())
tic
[detphi, phi] = rbrunforward(cfg);
toc

%# rbrunforward returned two outputs:
%# detphi - the sampled measurements (fluence, not diffuse reflectance!) at all src/det pairs
%# phi - the fluence at all nodes for all sources and detectors (column-dimension)

disp('== plot results ...')
figure;
colormap('jet')
plotmesh([cfg.node log10(phi(:,1))], cfg.elem, 'y > 25')  %# first column is from the src
%shading interp
title('forward solution from source')

figure
colormap('jet')
plotmesh([cfg.node log10(phi(:,2))], cfg.elem, 'y > 25')  %# 2nd column is from the detector
%shading interp
title('forward solution from detector')
