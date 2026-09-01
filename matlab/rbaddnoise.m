function newdata = rbaddnoise(data, snrshot, snrthermal, randseed)
%
% newdata=rbaddnoise(data, snrshot, snrthermal,randseed)
%
% Adding simulated shot-noise and thermal noise to the simulated data
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     data: noise-less data
%     snrshot: the desired SNR (in dB) of the shot-noise
%          (amplitude-proportional). the noise standard deviation follows
%          sigma(i)=sqrt(|data(i)|/max|data|)*max|data|*10^(-snrshot/20),
%          i.e. snrshot is the SNR of the strongest data point, and weaker
%          data points have a proportionally lower SNR, following
%              SNR(i) = snrshot + 10*log10(|data(i)|/max|data|)  (dB)
%          NOTE: this normalization was changed on 2026/09/01 - previously
%          the shot-noise was not referenced to max|data|, making the
%          requested SNR dependent on the units of the data. to convert an
%          snrshot value used with the old code, use
%              snrshot_new = snrshot_old + 10*log10(max|data|)
%     snrthermal: the desired SNR (in dB) of the thermal-noise (noise
%          floor), i.e. a data-independent sigma=max|data|*10^(-snrthermal/20)
%     randseed (optional): specify a seed for reproducible noise
%
%     for complex (frequency-domain) data, the noise is added as a
%     circularly-symmetric complex gaussian with sigma per quadrature, so
%     that the amplitude noise matches that of real-valued (CW) data at the
%     same snr, and the phase error follows sigma_phase(i)=sigma(i)/|data(i)|
%     (in radian) automatically. the total complex noise |newdata-data| is
%     therefore sqrt(2)*sigma.
%     if snrshot or snrthermal is given as a complex number, the imaginary
%     part specifies an additional independent phase jitter (instrument
%     phase noise) of 10^(-imag(snr)/20) radian, applied multiplicatively.
%     an imaginary part of 0 (i.e. a real-valued input) adds no jitter.
%
% output:
%     newdata: noise contaminated data
%
% license:
%     GPL version 3, see LICENSE_GPLv3.txt files for details
%
% -- this function is part of Redbird-m toolbox
%

if (nargin == 1)
    newdata = data;
    warning('no noise added');
    return
end

if (nargin < 4)
    randseed = 123456789;
end

if (exist('rng', 'file') == 2 || exist('rng', 'builtin') == 5)
    rng(randseed);
else
    randn('state', randseed);
end

if (nargin < 3)
    snrthermal = inf;
end

max_amp = max(abs(data(:)));

sigma_shot = sqrt(max_amp) * 10^(-real(snrshot) / 20);
sigma_thermal = max_amp * 10^(-real(snrthermal) / 20);

if (isreal(data))
    newdata = data + sqrt(abs(data)) .* randn(size(data)) * sigma_shot + ...
              randn(size(data)) * sigma_thermal;
else
    % circularly-symmetric complex gaussian, unit variance per quadrature,
    % so that the in-phase (amplitude) noise matches the real-data branch
    cnoise = @() complex(randn(size(data)), randn(size(data)));
    newdata = data + sqrt(abs(data)) .* cnoise() * sigma_shot + ...
              cnoise() * sigma_thermal;

    % optional extra phase jitter, requested via the imaginary part of the
    % snr inputs; independent contributions add in quadrature
    sigma_phase = 0;
    if (imag(snrshot) ~= 0)
        sigma_phase = sigma_phase + 10^(-imag(snrshot) / 10);
    end
    if (imag(snrthermal) ~= 0)
        sigma_phase = sigma_phase + 10^(-imag(snrthermal) / 10);
    end
    if (sigma_phase > 0)
        newdata = newdata .* exp(1i * randn(size(data)) * sqrt(sigma_phase));
    end
end
