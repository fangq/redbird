function [newJ, newy0, newphi] = rbmultispectral(sd, cfg, Jmua, y0, phi, params, rfcw, Jd, prop)
%
% [newJ, newy0, newphi]=rbmultispectral(Jmua, y0, phi, paramlist, Jd)
%
% Concatenate multi-special forward modeling data into a single linear
% system
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     Jmua: the multispectral Jacobian as a containers.Map object with each
%         wavelength as a key
%     y0: the multispectral measurement data
%     phi: the forward solution at multiple wavelengths
%     params: must be a struct as cfg.param
%     prop: must be a struct as cfg.prop
%
% output:
%     newJ: the contatenated Jacobian - dimension (Nw*Ns*Nd) x ((Nn*Np)+Nn*2), where
%          Nw - the number of wavelengths
%          Ns - the number of sources
%          Nd - the number of detectors
%          Nn - the number of nodes
%          Np - the number of parameter species
%          Nn*2 stores the scat-amplitude and scat-power if Jd is given
%     newy0: the contatenated measurement data vector
%     newphi: the contatenated forward simulation of measurements
%
% license:
%     GPL version 3, see LICENSE_GPLv3.txt files for details
%
% -- this function is part of Redbird-m toolbox
%

newJ = struct;
newy0 = [];
newphi = [];

ishelmholtz = isstruct(params) && (isfield(params, 'epsilon') || isfield(params, 'sigma'));

if (isa(Jmua, 'containers.Map') || (isstruct(Jmua) && isa(Jmua(1).J, 'containers.Map')))
    if (isa(Jmua, 'containers.Map'))
        wv = keys(Jmua);
    else
        wv = keys(Jmua(1).J);
        %         if (exist('Jd','var') && length(Jd) > 1)
        %             Jd = Jd(1).J;
        %         end
    end
    if (ishelmholtz)
        % MWT: chain Jmua -> Jeps, Jsigma per frequency
        % rbjac produces Jmua = -<E_s, E_r>_M (DOT sign convention).
        % For MWT (A has -k^2*M), Jk^2 = -Jmua. Then:
        %   Jeps   = (omega^2*mu0*eps0)  * Jk^2 = -(omega^2*mu0*eps0)  * Jmua
        %   Jsigma = (-1j*omega*mu0)     * Jk^2 = (1j*omega*mu0)       * Jmua
        eps0_mm = 8.854187817e-15;
        mu0_mm = 4 * pi * 1e-10;
        omegas = zeros(length(wv), 1);
        for ii = 1:length(wv)
            if isa(cfg.omega, 'containers.Map')
                omegas(ii) = cfg.omega(wv{ii});
            else
                omegas(ii) = cfg.omega;
            end
        end
        weight_eps = -(omegas.^2) * mu0_mm * eps0_mm;
        weight_sigma = 1j * omegas * mu0_mm;

        newJ = struct;
        for k = 1:rfcw
            if isa(Jmua, 'containers.Map')
                JJ = Jmua;
            else
                JJ = Jmua(k).J;
            end
            if (isfield(params, 'epsilon'))
                newJ(k).epsilon = rbmatflat(JJ, weight_eps);
            end
            if (isfield(params, 'sigma'))
                newJ(k).sigma = rbmatflat(JJ, weight_sigma);
            end
        end
        if (length(newJ) == 1)
            newJ = newJ(1);
        end
        % handle y0/phi flattening (shared with DOT path below)
    else
        paramlist = fieldnames(params);
        has_norm_scat = (length(intersect(paramlist, {'scatamp500', 'scatpow500'})) == 2);
        has_legacy_scat = (length(intersect(paramlist, {'scatamp', 'scatpow'})) == 2);
        if (nargin > 7 && (has_norm_scat || has_legacy_scat))
            dcoeff = containers.Map();
            for i = 1:length(wv)
                dtemp = prop(wv{i});
                if ((isa(Jd, 'containers.Map') && size(dtemp, 1) < size(Jd(wv{i}), 2)) || ((isstruct(Jd) && size(dtemp, 1) < size(Jd(1).J(wv{i}), 2))))
                    dtemp = dtemp(2:end, :);
                end
                dcoeff(wv{i}) = 1 ./ (3 .* (dtemp(:, 1) + dtemp(:, 2)))';
            end
            if (has_norm_scat)
                Jscat = rbjacscat(Jd, dcoeff, params.scatpow500, wv, 500, '500');
            else
                Jscat = rbjacscat(Jd, dcoeff, params.scatpow, wv);
            end
        end
        chromophores = intersect(paramlist, {'hbo', 'hbr', 'water', 'lipids', 'aa3'});

        newJ = rbjacchrome(Jmua, chromophores);
        if (exist('Jscat', 'var') && isstruct(Jscat))
            allkeys = fieldnames(Jscat);
            for i = 1:length(allkeys)
                if (isstruct(newJ) && length(newJ) > 1)
                    newJ(1).(allkeys{i}) = Jscat(1).(allkeys{i});
                    newJ(2).(allkeys{i}) = Jscat(2).(allkeys{i});
                else
                    newJ.(allkeys{i}) = Jscat.(allkeys{i});
                end
            end
            clear Jscat;
        end
    end  % end of if (ishelmholtz) ... else (DOT chromophore path)
else
    newJ.mua = Jmua;
    if (nargin > 7 && ~isa(Jd, 'containers.Map'))
        newJ.dcoeff = Jd;
    end
end

if (nargin > 7 && ~isa(Jd, 'containers.Map') && ~isstruct(Jd))
    newJ.dcoeff = Jd;
end

if (~isa(y0, 'containers.Map') && ~isfield(y0, 'detphi'))
    newy0 = y0;
else
    if isstruct(y0)
        wv = keys(y0(1).detphi);
        %         rfcw = length(phi);
    else
        wv = keys(y0);
        %         rfcw = 1;
        temp = struct('y0', y0);
        clear y0;
        y0(rfcw).detphi = temp.y0;
        clear temp;
    end
    newy0 = struct('detphi', cell(1, max(rfcw)));
    for i = 1:length(wv)
        sdwv = sd(wv{i});
        if (size(sdwv, 2) == 3)
            sdwv(:, 4) = rfcw;
        end
        for j = rfcw
            tempphi = reshape(y0(j).detphi(wv{i}), [], 1);
            sdtemp = sdwv(sdwv(:, 4) == j | sdwv(:, 4) == 3, :);
            tempphi = tempphi(find(sdtemp(:, 3)));
            newy0(j).detphi = [newy0(j).detphi; tempphi];
            clear tempphi;
        end
    end

    if (length(rfcw) == 1)
        newy0 = newy0(rfcw).detphi;
    end
end

if (~isa(phi, 'containers.Map') && ~isfield(phi, 'detphi'))
    newphi = phi;
else
    if isstruct(phi)
        wv = keys(phi(1).detphi);
        %         rfcw = length(phi);
    else
        wv = keys(phi);
        %         rfcw = 1;
        temp = struct('phi', phi);
        clear phi;
        phi(rfcw).detphi = temp.phi;
        clear temp;
    end
    newphi = struct('detphi', cell(1, max(rfcw)));
    for i = 1:length(wv)
        sdwv = sd(wv{i});
        if (size(sdwv, 2) == 3)
            sdwv(:, 4) = rfcw;
        end
        for j = rfcw
            tempphi = reshape(phi(j).detphi(wv{i}), [], 1);
            sdtemp = sdwv(sdwv(:, 4) == j | sdwv(:, 4) == 3, :);
            tempphi = tempphi(find(sdtemp(:, 3)));
            newphi(j).detphi = [newphi(j).detphi; tempphi];
            clear tempphi sdtemp;
        end
    end

    if (length(rfcw) == 1)
        newphi = newphi(rfcw).detphi;
    end
end
