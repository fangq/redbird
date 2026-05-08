function [Amat, deldotdel] = rbfemlhs(cfg, deldotdel, wavelength, mode)
%
% [Amat,deldotdel]=rbfemlhs(cfg)
%   or
% [Amat,deldotdel]=rbfemlhs(cfg, wavelength)
% [Amat,deldotdel]=rbfemlhs(cfg, deldotdel, wavelength)
%
% create the FEM stiffness matrix (left-hand-side) for solving the
% diffusion equation (DOT) or the scalar Helmholtz equation (MWT) at a
% given wavelength/frequency. MWT mode is engaged when cfg.bulk.epsilon
% or cfg.bulk.sigma is present (defining the bulk medium for the
% first-order Bayliss-Turkel radiation boundary condition).
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     cfg: the initial simulation data structure
%     deldotdel (optional): precomputed operator on the mesh (del_phi dot del_phi)
%         where del represents the gradient; see help rbdeldotdel
%     wavelength (optional): a string or number denoting the wavelength
%
% output:
%     Amat: the left-hand-side matrix of the FEM equation - a sparse matrix
%          of dimension Nn x Nn, where Nn is the number of nodes of the
%          forward mesh
%     deldotdel: if the 2nd input is not given, this function compute
%          deldotdel and return as the 2nd output
%
% license:
%     GPL version 3, see LICENSE_GPLv3.txt files for details
%
% -- this function is part of Redbird-m toolbox
%

nn = size(cfg.node, 1);
ne = size(cfg.elem, 1);

R_C0 = (1 ./ 299792458000.);

% MWT (Helmholtz) detection: cfg.bulk.epsilon or cfg.bulk.sigma defines the
% bulk medium for the radiation boundary condition.
ishelmholtz = isfield(cfg, 'bulk') && (isfield(cfg.bulk, 'epsilon') || isfield(cfg.bulk, 'sigma'));

% cfg.prop is updated from cfg.param and contains the updated material props.
% DOT: cfg.prop columns are [mua mus g n] (1/mm, 1/mm, scalar, scalar)
% MWT: cfg.prop columns are [eps_r sigma mu0 n] (scalar, S/mm, H/mm, scalar)
% if cfg.param is node/elem based, cfg.prop is updated accordingly.
% if cfg.param is segmentation based, cfg.prop has the same format as mcx's
% prop, where the first row is label 0, and total length is Nseg+1

prop = cfg.prop;
if (~ishelmholtz)
    cfgreff = cfg.reff;
end
if exist('mode', 'var')
    if mode == 1
        omega = cfg.omega;
    elseif mode == 2
        omega = 0;
    end
else
    omega = cfg.omega;
end

if (nargin == 2 && numel(deldotdel) == 1)
    wavelength = deldotdel;
end

if (isa(cfg.prop, 'containers.Map')) % if multiple wavelengths, take current
    if (nargin < 3)
        error('you must specify wavelength');
    end
    if (~ischar(wavelength))
        wavelength = sprintf('%g', wavelength);
    end
    prop = cfg.prop(wavelength);
    if (~ishelmholtz)
        cfgreff = cfg.reff(wavelength);
    end
    if (isa(omega, 'containers.Map'))
        omega = omega(wavelength);
    end
end

% if deldotdel is provided, call native code; otherwise, call mex

if (nargin >= 2 && numel(deldotdel) > 1)
    % extract material properties for the FEM volume assembly. After this
    % block, (avol, breal, bimag) define the per-element or per-node
    % coefficients in:  A_e = avol*<grad phi_i, grad phi_j>_e
    %                       + (breal + 1j*bimag)*<phi_i, phi_j>_e
    if (ishelmholtz)
        eps0_mm = 8.854187817e-15;  % F/mm
        if (size(prop, 1) == nn || size(prop, 1) == ne)
            eps_r = prop(:, 1);
            sigma = prop(:, 2);
            permea = prop(:, 3);
        elseif (size(prop, 1) < min([nn ne])) % segmentation based prop list
            eps_r = prop(cfg.seg + 1, 1);
            sigma = prop(cfg.seg + 1, 2);
            permea = prop(cfg.seg + 1, 3);
        end
        % k^2 (1/mm^2) = w^2*mu*eps0*eps_r - j*w*mu*sigma. Mass coefficient
        % in the weak form is -k^2.
        avol = ones(size(eps_r));
        breal = -(omega.^2) .* permea .* eps0_mm .* eps_r;
        bimag = omega .* permea .* sigma;
    else
        if (size(prop, 1) == nn || size(prop, 1) == ne)
            mua = prop(:, 1);
            if (size(prop, 2) < 3)
                musp = prop(:, 2);
            else
                musp = prop(:, 2) .* (1 - prop(:, 3));
            end
        elseif (size(prop, 1) < min([nn ne])) % use segmentation based prop list
            mua = prop(cfg.seg + 1, 1);
            if (size(prop, 2) < 3)
                musp = prop(cfg.seg + 1, 2); % assume g is 0
            else
                musp = prop(cfg.seg + 1, 2) .* (1 - prop(cfg.seg + 1, 3));
            end
        end
        if (isfield(cfg, 'bulk') && isfield(cfg.bulk, 'n'))
            nref = cfg.bulk.n;
        elseif (isfield(cfg, 'seg') && size(prop, 1) < min([nn, ne]))
            nref = rbgetbulk(cfg);
            if (isa(nref, 'containers.Map'))
                nref = nref(wavelength);
            end
            nref = nref(4);
        else
            nref = prop(:, 4);
        end
        if (length(nref) == 1)
            nref = nref * ones(size(mua));
        end
        avol = 1 ./ (3 .* (mua + musp));
        breal = mua;
        bimag = omega .* R_C0 .* nref;
    end

    edges = sort(meshedge(cfg.elem), 2);

    % volume assembly: dispatches on per-element vs per-node properties
    if (length(breal) == size(cfg.elem, 1))  % element based property
        Aoffd = deldotdel(:, [2:4, 6:7, 9]) .* repmat(avol(:), 1, 6) + repmat(0.05 .* breal(:) .* cfg.evol(:), 1, 6);
        Adiag = deldotdel(:, [1, 5, 8, 10]) .* repmat(avol(:), 1, 4) + repmat(0.10 .* breal(:) .* cfg.evol(:), 1, 4);
        if (any(bimag(:) ~= 0))
            Aoffd = complex(Aoffd, repmat(0.05 .* bimag(:) .* cfg.evol(:), 1, 6));
            Adiag = complex(Adiag, repmat(0.10 .* bimag(:) .* cfg.evol(:), 1, 4));
        end
    else  % node based properties
        w1 = (1 / 120) * [2 2 1 1; 2 1 2 1; 2 1 1 2; 1 2 2 1; 1 2 1 2; 1 1 2 2]';
        w2 = (1 / 60) * (diag([2 2 2 2]) + 1);
        breal_e = reshape(breal(cfg.elem), size(cfg.elem));
        bimag_e = reshape(bimag(cfg.elem), size(cfg.elem));
        avol_e = mean(reshape(avol(cfg.elem), size(cfg.elem)), 2);
        Aoffd = deldotdel(:, [2:4, 6:7, 9]) .* repmat(avol_e, 1, 6) + (breal_e * w1) .* repmat(cfg.evol(:), 1, 6);
        Adiag = deldotdel(:, [1, 5, 8, 10]) .* repmat(avol_e, 1, 4) + (breal_e * w2) .* repmat(cfg.evol(:), 1, 4);
        if (any(bimag(:) ~= 0))
            Aoffd = complex(Aoffd, (bimag_e * w1) .* repmat(cfg.evol(:), 1, 6));
            Adiag = complex(Adiag, (bimag_e * w2) .* repmat(cfg.evol(:), 1, 4));
        end
    end

    % boundary condition: Robin (DOT) or first-order Bayliss-Turkel (MWT)
    edgebc = sort(meshedge(cfg.face), 2);
    if (ishelmholtz)
        bk = rbgetbulk(cfg);
        if (isa(bk, 'containers.Map'))
            bk = bk(wavelength);
        end
        k2bg = (omega.^2) * bk(3) * eps0_mm * bk(1) - 1j * omega * bk(3) * bk(2);
        kbg = sqrt(k2bg);
        rvec = cfg.facecenter - repmat(cfg.rbcorigin, size(cfg.facecenter, 1), 1);
        rdotn = sum(rvec .* cfg.facenormal, 2) ./ cfg.facer;
        bccoef = (1j * kbg - 1 ./ (2 * cfg.facer)) .* rdotn;
        Adiagbc = (cfg.area(:) / 6) .* bccoef;
    else
        Reff = cfgreff;
        Adiagbc = cfg.area(:) * ((1 - Reff) / (12 * (1 + Reff)));
    end
    Adiagbc = repmat(Adiagbc, 1, 3);
    Aoffdbc = Adiagbc * 0.5;

    Amat = sparse([edges(:, 1); edges(:, 2); cfg.elem(:); edgebc(:, 1); edgebc(:, 2); cfg.face(:)], ...
                  [edges(:, 2); edges(:, 1); cfg.elem(:); edgebc(:, 2); edgebc(:, 1); cfg.face(:)], ...
                  [Aoffd(:); Aoffd(:); Adiag(:); Aoffdbc(:); Aoffdbc(:); Adiagbc(:)]);
else
    if (size(cfg.elem, 2) > 4)
        cfg.elem(:, 5:end) = [];
    end
    cfg.prop = prop; % use property of the current wavelength
    [Adiag, Aoffd, deldotdel] = rbfemmatrix(cfg);
    Amat = sparse([cfg.rows, cfg.cols, (1:nn)], [cfg.cols, cfg.rows, (1:nn)], [Aoffd, Aoffd, Adiag], nn, nn);
    deldotdel = deldotdel';
end
