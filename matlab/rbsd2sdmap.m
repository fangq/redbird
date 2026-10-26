function [sdmap, sdidx] = rbsd2sdmap(cfg, sd, md)
%
% [sdmap, sdidx] = rbsd2sdmap(cfg, sd)
%   or
% [sdmap, sdidx] = rbsd2sdmap(cfg, sd, md)
%
% Convert a Redbird source-detector table into the cfg.sdmap input of
% mcxlab/mmclab, which then compute the adjoint Jacobian of these channels only
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     cfg: the redbird data structure, defining srcpos, detpos and optionally widesrc
%     sd: the source-detector table of one wavelength, as returned by rbsdmap:
%         [src det active mode]; det counts after all sources (srcnum+widesrcnum+d)
%     md: the measurement mode of the Jacobian, 1: RF/CW (default), 2: CW;
%         rows with sd(:,4) equal to md or 3 are kept, as rbjac does
%
% output:
%     sdmap: an Nch x 2 array of 1-based [srcid detid] rows, detid counting the
%         rows of cfg.detpos, for the active rows of sd in their order; the rows
%         of the Monte Carlo Jacobian then match those of rbjac
%     sdidx: the index of each channel in detphi(:), where detphi is the
%         Ndet x Nsrc detector reading matrix returned by rbrunforward
%
% license:
%     GPL version 3, see LICENSE_GPLv3.txt files for details
%
% -- this function is part of Redbird toolbox
%

if (nargin < 3)
    md = 1;
end

srcnum = size(cfg.srcpos, 1);
detnum = size(cfg.detpos, 1);
detoffset = srcnum;
if (isfield(cfg, 'widesrc'))
    detoffset = detoffset + size(cfg.widesrc, 1);
end

if (size(sd, 2) < 4)
    sd(:, 4) = md;
end

sd = sd(sd(:, 3) == 1 & (sd(:, 4) == md | sd(:, 4) == 3), :);

if (any(sd(:, 1) > srcnum | sd(:, 2) <= detoffset | sd(:, 2) > detoffset + detnum))
    error('rbsd2sdmap: wide-field sources or detectors are not supported by the Monte Carlo Jacobian');
end

sdmap = [sd(:, 1), sd(:, 2) - detoffset];
sdidx = (sdmap(:, 1) - 1) * detnum + sdmap(:, 2);
