function [rhs, loc, bary, optode] = rbfemrhs(cfg, sd, wv, md)
%
% newcfg=rbmeshprep(cfg)
%
% Create the right-hand-sides for the FEM system equation, here we solve
% forward systems for both source and detector locations in order to use
% the adjoint method to create Jacobians.
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     cfg: the initial simulation data structure
%
% output:
%     rhs: the right-hand-side of the FEM equation; if multiple sources are
%          used, the rhs is a matrix of dimension Nn x (Ns+Nd), where Nn is
%          the total number of nodes, and Ns the number of sources and
%          Nd is the total number of detectors
%     loc: the indices of the forward mesh element that encloses each
%          source or detector; Nan means the source is outside of the mesh
%          or a wide-field source/detector
%     bary: the barycentric coordinates of the source/detector if it is
%          enclosed by a tetrahedral element; Nan if outside of the mesh
%     optode: the full source+detector position list returned by
%          rbgetoptodes.m
%
% license:
%     GPL version 3, see LICENSE_GPLv3.txt files for details
%
% -- this function is part of Redbird toolbox
%

if (nargin > 2)
    sd = sd(wv);
    if (nargin > 3 && size(sd, 2) > 3)
        sd = sd((sd(:, 4) == md | sd(:, 4) == 3), :);
        if (isfield(cfg, 'rhsweight') && isstruct(cfg.rhsweight))
            rhsweight = cfg.rhsweight(md).weight;
        end
    end
    [optsrc, optdet, widesrc, widedet] = rbgetoptodes(cfg, wv);
else
    [optsrc, optdet, widesrc, widedet] = rbgetoptodes(cfg);
end

if exist('sd', 'var')
    if (size(sd, 2) > 3)
        tempsd = sd(sd(:, 3) == 1 & ismember(sd(:, 4), [1 2 3]), :);
    else
        tempsd = sd(sd(:, 3) == 1, :);
    end
    activeOpt = unique(tempsd(:, 1:2));
else
    activeOpt = 1:size([optsrc; optdet], 1) + size([widesrc; widedet], 1);
end

if ((size(optsrc, 1) < 1 && size(widesrc, 1) < 1) || (size(optdet, 1) < 1 && size(widedet, 1) < 1))
    error('you must provide at least one source or detector');
end

loc = [];
bary = [];

rhs = sparse(size(cfg.node, 1), size([widesrc; widedet], 1) + size([optsrc; optdet], 1));
if ~isempty(optsrc)
    [locsrc, barysrc] = tsearchn(cfg.node, cfg.elem, optsrc);
    rhsoptsrc = point2rhs(optsrc, locsrc, barysrc, cfg);
else
    locsrc = [];
    barysrc = [];
    rhsoptsrc = [];
end
if ~isempty(optdet)
    [locdet, barydet] = tsearchn(cfg.node, cfg.elem, optdet);
    rhsoptdet = point2rhs(optdet, locdet, barydet, cfg);
else
    locdet = [];
    barydet = [];
    rhsoptdet = [];
end

loc = [locsrc; nan * ones(size(widesrc, 1), 1); locdet; nan * ones(size(widedet, 1), 1)];
bary = [barysrc; nan * ones(size(widesrc, 1), 4); barydet; nan * ones(size(widesrc, 1), 4)];

rhs = [rhsoptsrc widesrc' rhsoptdet widedet'];

inactive = setdiff(1:size(rhs, 2), activeOpt);

%% Test rhs weight
if (isfield(cfg, 'rhsweight') && ~exist('rhsweight', 'var'))
    rhs = rhs .* cfg.rhsweight.';
elseif exist('rhsweight', 'var')
    rhs = rhs .* rhsweight.';
end
%%

rhs(:, inactive) = 0;

function optoderhs = point2rhs(srcs, locidx, baryc, cfg)
if (size(srcs, 1) ~= size(baryc, 1))
    error('barycentric indices should equal number of sources or detectors');
end

optoderhs = zeros(size(cfg.node, 1), size(srcs, 1));
for ii = 1:size(srcs, 1)
    if ~isnan(locidx(ii))
        optoderhs(cfg.elem(locidx(ii), :), ii) = baryc(ii, :)';
    end
end
