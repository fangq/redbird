clear cfg

disp('== mesh generation (iso2mesh) ...')
%# Create a 3 layer domain
maxvol = 4; %# max tet element volume

tic
[no1, fc1, regionseeds] = latticegrid([0,60], [0, 50], [0, 5, 10, 30]);%# create a 3-layer z-lattice
[cfg.node, cfg.elem] = s2m(no1,fc1,1,maxvol, 'tetgen', regionseeds);   %# generate tetrahedral mesh, label increases in z-axis

cfg.seg = cfg.elem(:,5);    %# cfg.seg is similar to cfg.elemprop in mmclab, but also supports node-based labels
cfg.elem(:,5)=[];
toc

plotmesh(cfg.node, [cfg.elem, cfg.seg], 'y>20')

disp('== define settings in cfg ...')

%# Creating forward simulation data structure cfg
%# properties:[mua(1/mm),mus(1/mm),g, n]  % if both mus/g are given, mus'=mus/(1-g) will be used for diffusion, usually set g to 0
cfg.prop  =  [0,        0,        1, 1           %# cfg.prop is the same as mcx/mmc, row-1 is for label 0, row-2 for label 1 etc
              0.006,    0.8,      0, 1.37        %# label 1 (0 < z < 5)
              0.02,     0.4,      0, 1.37        %# label 2 (5 < z < 10)
              0.002,    1,        0, 1.37];      %# label 2 (10 < z < 30)

%# Default redbird source is a pencil beam, same as mcx/mmc
[xi, yi] = meshgrid(5:5:50, 5:5:40);
cfg.srcpos = [xi(:), yi(:), ones(length(xi(:)), 1)];    %# redbird srcpos can have multiple rows
cfg.srcdir = [0, 0, 1];          %# srcdir determines how sources are sunken into the mesh

%# Redbird detector positions are point-like, directly sampling the output fluence; it is different from the disk-like shape in mcx/mmc
[xi, yi] = meshgrid(7.5:5:50, 7.5:5:40);
cfg.detpos = [xi(:), yi(:), ones(length(xi(:)), 1)];  %# redbird detpos can have multiple rows
cfg.detdir = [0, 0, 1];          %# redbird automatically computes the adjoint solutions, treating detectors as sources

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
plotmesh([cfg.node log10(abs(phi(:,1)))], cfg.elem, 'y > 25')  %# first column is from the src
%shading interp
title('forward solution from source')

figure
colormap('jet')
plotmesh([cfg.node log10(abs(phi(:,2)))], cfg.elem, 'y > 25')  %# 2nd column is from the detector
%shading interp
title('forward solution from detector')
