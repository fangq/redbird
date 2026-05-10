function [newcfg, sd] = rbmeshprep(cfg)
%
% newcfg=rbmeshprep(cfg)
%
% Compute all missing fields from the cfg input sturcture to get
% ready for forward and inverse modeling
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     cfg: the initial simulation data structure
%
% output:
%     newcfg: the updated simulation data structure after adding all missing fields
%
% license:
%     GPL version 3, see LICENSE_GPLv3.txt files for details
%
% -- this function is part of Redbird-m toolbox
%

if (~isfield(cfg, 'node') || ~isfield(cfg, 'elem'))
    error('cfg.node or cfg.elem is missing');
end

cfg.elem(:, 1:4) = meshreorient(cfg.node(:, 1:3), cfg.elem(:, 1:4));

if ((~isfield(cfg, 'seg') || isempty(cfg.seg)) && size(cfg.elem, 2) > 4)
    cfg.seg = cfg.elem(:, 5);
    cfg.elem(:, 5) = [];
end
if (~isfield(cfg, 'isreoriented') || isempty(cfg.isreoriented) || cfg.isreoriented == 0)
    cfg.elem = meshreorient(cfg.node, cfg.elem(:, 1:4));
    cfg.isreoriented = 1;
end
if (~isfield(cfg, 'face') || isempty(cfg.face))
    cfg.face = volface(cfg.elem);
end
if (~isfield(cfg, 'area') || isempty(cfg.area))
    cfg.area = elemvolume(cfg.node, cfg.face);
end
if (~isfield(cfg, 'evol') || isempty(cfg.evol))
    cfg.evol = elemvolume(cfg.node, cfg.elem);
end
if (~isfield(cfg, 'nvol') || isempty(cfg.nvol))
    cfg.nvol = nodevolume(cfg.node, cfg.elem, cfg.evol);
end
if (find(cfg.evol == 0))
    fprintf(1, ['degenerated elements are detected: [' sprintf('%d ', find(cfg.evol == 0)) ']\n']);
    error(['input mesh can not contain degenerated elements, ' ...
           'please double check your input mesh; if you use a ' ...
           'widefield source, please rerun mmcsrcdomain and setting ' ...
           '''Expansion'' option to a larger value (default is 1)']);
end
if (~isfield(cfg, 'srcpos'))
    error('cfg.srcpos field is missing');
end
src_is_line = isfield(cfg, 'srctype') && strcmp(cfg.srctype, 'line');
src_is_ray  = isfield(cfg, 'srctype') && strcmp(cfg.srctype, 'ray');
if (~isfield(cfg, 'srcdir') && ~src_is_line)
    error('cfg.srcdir field is missing');
end
if (isfield(cfg, 'prop') && isfield(cfg, 'param') && ...
    isa(cfg.prop, 'containers.Map'))
    wv = cfg.prop.keys;
    if (~isempty(wv))
        cfg.prop = rbupdateprop(cfg);
    end
end

ishelmholtz = isfield(cfg, 'bulk') && (isfield(cfg.bulk, 'epsilon') || isfield(cfg.bulk, 'sigma'));

% face-adjacency table: 4 face-neighbors per tet (0 = boundary face).
% used by line-source tracing for both DOT and MWT.
if (~isfield(cfg, 'facenb') || isempty(cfg.facenb))
    cfg.facenb = faceneighbors(cfg.elem(:, 1:4));
end

if (~ishelmholtz)
    % compute R_eff - effective reflection coeff, and musp0 - background mus'
    if (~isfield(cfg, 'reff') || isempty(cfg.reff))
        bkprop = rbgetbulk(cfg);
        if (isa(bkprop, 'containers.Map'))
            cfg.reff = containers.Map();
            cfg.musp0 = containers.Map();
            for waveid = bkprop.keys
                wv = waveid{1};
                prop = bkprop(wv);
                cfg.reff(wv) = rbgetreff(prop(4), 1);
                cfg.musp0(wv) = prop(2) * (1 - prop(3));
            end
        else
            cfg.reff = rbgetreff(bkprop(4), 1);
            cfg.musp0 = bkprop(2) * (1 - bkprop(3));
        end
    end
else
    % MWT: precompute Bayliss-Turkel RBC geometry on cfg.face
    if (~isfield(cfg, 'facecenter') || isempty(cfg.facecenter))
        cfg.facecenter = (cfg.node(cfg.face(:, 1), 1:3) + ...
                          cfg.node(cfg.face(:, 2), 1:3) + ...
                          cfg.node(cfg.face(:, 3), 1:3)) / 3;
    end
    if (~isfield(cfg, 'facenormal') || isempty(cfg.facenormal))
        ab = cfg.node(cfg.face(:, 2), 1:3) - cfg.node(cfg.face(:, 1), 1:3);
        ac = cfg.node(cfg.face(:, 3), 1:3) - cfg.node(cfg.face(:, 1), 1:3);
        nrm = cross(ab, ac, 2);
        nlen = sqrt(sum(nrm .* nrm, 2));
        cfg.facenormal = nrm ./ repmat(nlen, 1, 3);
    end
    if (~isfield(cfg, 'rbcorigin') || isempty(cfg.rbcorigin))
        if (isfield(cfg, 'srcpos') && isfield(cfg, 'detpos'))
            optodes = [cfg.srcpos(:, 1:3); cfg.detpos(:, 1:3)];
            if (isfield(cfg, 'srctype') && strcmp(cfg.srctype, 'line') && isfield(cfg, 'srcparam1') && numel(cfg.srcparam1) >= 3)
                optodes = [optodes; cfg.srcpos(:, 1:3) + repmat(cfg.srcparam1(1:3), size(cfg.srcpos, 1), 1)];
            end
            if (isfield(cfg, 'dettype') && strcmp(cfg.dettype, 'line') && isfield(cfg, 'detparam1') && numel(cfg.detparam1) >= 3)
                optodes = [optodes; cfg.detpos(:, 1:3) + repmat(cfg.detparam1(1:3), size(cfg.detpos, 1), 1)];
            end
            cfg.rbcorigin = mean(optodes, 1);
        else
            cfg.rbcorigin = mean(cfg.node(:, 1:3), 1);
        end
    end
    if (~isfield(cfg, 'facer') || isempty(cfg.facer))
        rvec = cfg.facecenter - repmat(cfg.rbcorigin, size(cfg.facecenter, 1), 1);
        cfg.facer = sqrt(sum(rvec .* rvec, 2));
    end
end
det_is_line = isfield(cfg, 'detpos') && isfield(cfg, 'dettype') && strcmp(cfg.dettype, 'line');
det_is_ray  = isfield(cfg, 'detpos') && isfield(cfg, 'dettype') && strcmp(cfg.dettype, 'ray');
if (((isfield(cfg, 'srctype') && ~ismember(cfg.srctype, {'pencil', 'isotropic'})) || isfield(cfg, 'widesrcid') || src_is_line || src_is_ray) && ~isfield(cfg, 'widesrc'))
    cfg.srcpos0 = cfg.srcpos;
    cfg = rbsrc2bc(cfg);
end
if (((isfield(cfg, 'dettype') && ~ismember(cfg.dettype, {'pencil', 'isotropic'})) || isfield(cfg, 'widedetid') || det_is_line || det_is_ray) && ~isfield(cfg, 'widedet'))
    cfg.detpos0 = cfg.detpos;
    cfg = rbsrc2bc(cfg, 1);
end
if (~isfield(cfg, 'cols') || isempty(cfg.cols))
    [cfg.rows, cfg.cols, cfg.idxcount] = rbfemnz(cfg.elem, size(cfg.node, 1));
end

if (~isfield(cfg, 'idxsum') || isempty(cfg.idxsum))
    cfg.idxsum = cumsum(cfg.idxcount);
end

if (~isfield(cfg, 'deldotdel') || isempty(cfg.deldotdel))
    cfg.deldotdel = rbdeldotdel(cfg);
end
if (~isfield(cfg, 'omega') || isempty(cfg.omega))
    cfg.omega = 0;
end

% time-domain DOT (Crank-Nicolson) is engaged when cfg.tstart, cfg.tstep,
% and cfg.tend are all defined. Validate and reject mixed configurations.
istd = isfield(cfg, 'tstart') && ~isempty(cfg.tstart) && ...
       isfield(cfg, 'tstep') && ~isempty(cfg.tstep) && ...
       isfield(cfg, 'tend') && ~isempty(cfg.tend);
if (istd)
    if (cfg.tend <= cfg.tstart)
        error('cfg.tend must be greater than cfg.tstart for time-domain DOT');
    end
    if (cfg.tstep <= 0)
        error('cfg.tstep must be positive for time-domain DOT');
    end
    if (any(cfg.omega(:) ~= 0))
        error('time-domain DOT (cfg.tstart/tstep/tend) cannot be used with frequency-domain modulation (cfg.omega>0); pick one');
    end
    if (ishelmholtz)
        error('time-domain Crank-Nicolson is only valid for the parabolic diffusion equation; remove cfg.bulk.epsilon/sigma to use TD');
    end
end
newcfg = cfg;
sd = rbsdmap(cfg);
