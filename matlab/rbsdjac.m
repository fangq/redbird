function [Asd, srcch, detch] = rbsdjac(sd, wvkeys, detphiflat, csd, sdindex, reform, sdmode, rfcw)
%
% [Asd, srcch, detch] = rbsdjac(sd, wvkeys, detphiflat, csd, sdindex, reform, sdmode, rfcw)
%
% Build the real-valued Jacobian block for simultaneous source-detector
% (SD) coupling coefficient estimation. The measurement model is
%     phi_hat(i,j) = csd(srcch(i)) * csd(detch(j)) * phi(i,j)
% so the complex Jacobian entries are J(:,ch) = model/csd(ch) for each of
% the two roles (source and detector) of every measurement. See the
% "simultaneous source-detector coupling coefficient estimation" section
% of the Redbird manual and Fang 2004 (Redbird90 rb2_inversion.f90).
%
% input:
%     sd: the source-detector mapping table (containers.Map or array);
%         rows must be ordered consistently with the flattened data
%     wvkeys: cell array of wavelength keys defining the row order of the
%         flattened data (use keys(detphi0)); {''} for a plain sd table
%     detphiflat: the flattened complex model prediction (already scaled
%         by the current coupling coefficients), one entry per active
%         measurement, concatenated over wavelengths
%     csd: current complex coupling coefficient vector (one per channel)
%     sdindex: maps each optode (sources first, then detectors, using the
%         same global indices as the sd table columns 1-2) to a coupling
%         channel; tie a source and detector to the same channel when they
%         are the same physical optode (e.g. MWT antennas)
%     reform: 'real'/'reim' for the Cauchy-Riemann expanded complex form,
%         or 'logphase' for the approximate log-magnitude/phase form
%         (J_logamp=1/|csd|, J_phase=1, see manual Eq. J5/J6)
%     sdmode: 'y' (default) estimate both roles, 's' source-only,
%         'd' detector-only
%     rfcw: (optional) source modality flag as in rbrunrecon, default 1
%
% output:
%     Asd: the real Jacobian block, size 2*Nmeas x 2*Nchannel; the first
%         Nchannel columns are Re(dc) (or d|c| for logphase), the last
%         Nchannel columns Im(dc) (or the phase of c)
%     srcch, detch: per-measurement-row coupling channel indices
%
% license:
%     GPL version 3, see LICENSE_GPLv3.txt files for details
%
% -- this function is part of Redbird toolbox
%

if (nargin < 8)
    rfcw = 1;
end
if (nargin < 7 || isempty(sdmode))
    sdmode = 'y';
end

nch = max(sdindex);
srcch = [];
detch = [];

for i = 1:length(wvkeys)
    if (isa(sd, 'containers.Map'))
        sdwv = sd(wvkeys{i});
    else
        sdwv = sd;
    end
    if (size(sdwv, 2) > 3)
        sdwv = sdwv((sdwv(:, 4) == rfcw | sdwv(:, 4) == 3), :);
    end
    kept = (sdwv(:, 3) == 1);
    srcch = [srcch; sdindex(sdwv(kept, 1))'];
    detch = [detch; sdindex(sdwv(kept, 2))'];
end
srcch = srcch(:);
detch = detch(:);

nrows = length(detphiflat);
if (length(srcch) ~= nrows)
    error('sd table rows (%d) do not match the flattened data length (%d)', length(srcch), nrows);
end

fitsrc = (sdmode == 'y' || sdmode == 's');
fitdet = (sdmode == 'y' || sdmode == 'd');
rows = (1:nrows)';

if (strcmp(reform, 'logphase'))
    % approximate log-amplitude/phase SD Jacobian: d(log|phi_hat|)/d|c| = 1/|c|,
    % d(angle(phi_hat))/d(angle(c)) = 1; cross-terms are zero
    Alog = zeros(nrows, nch);
    Apha = zeros(nrows, nch);
    if (fitsrc)
        Alog(sub2ind([nrows, nch], rows, srcch)) = 1 ./ abs(csd(srcch));
        Apha(sub2ind([nrows, nch], rows, srcch)) = 1;
    end
    if (fitdet)
        Alog(sub2ind([nrows, nch], rows, detch)) = Alog(sub2ind([nrows, nch], rows, detch)) + 1 ./ abs(csd(detch));
        Apha(sub2ind([nrows, nch], rows, detch)) = Apha(sub2ind([nrows, nch], rows, detch)) + 1;
    end
    Asd = [Alog, zeros(nrows, nch); zeros(nrows, nch), Apha];
else
    % complex form: J(:,ch)=model/csd(ch), expanded as [Re -Im; Im Re] so the
    % unknowns are [Re(dc); Im(dc)]
    Jsd = zeros(nrows, nch);
    if (fitsrc)
        Jsd(sub2ind([nrows, nch], rows, srcch)) = detphiflat(:) ./ csd(srcch);
    end
    if (fitdet)
        Jsd(sub2ind([nrows, nch], rows, detch)) = Jsd(sub2ind([nrows, nch], rows, detch)) + detphiflat(:) ./ csd(detch);
    end
    Asd = [real(Jsd), -imag(Jsd); imag(Jsd), real(Jsd)];
end
