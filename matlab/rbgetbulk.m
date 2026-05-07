function bkprop = rbgetbulk(cfg)
%
% bkprop=rbgetbulk(cfg)
%
% Return the optical properties of the "bulk" medium, which is considered
% the medium on the outer-most layer and is interfaced directly with air
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     cfg: the forward simulation data structure
%
% output:
%     bkprop: the optical property quadruplet in the order of
%             [mua(1/mm), mus(1/mm), g, n]
%         if single wavelength, and a containers.Map object made of the
%         above quadruplet for each wavelength.
%
% license:
%     GPL version 3, see LICENSE_GPLv3.txt files for details
%
% -- this function is part of Redbird-m toolbox
%

bkprop = [0 0 0 1.37];

if (~isfield(cfg, 'bulk'))
    if (isfield(cfg, 'prop') && ~isempty(cfg.prop))
        nn = size(cfg.node, 1);
        ne = size(cfg.elem, 1);
        prop = containers.Map();
        bkprop = containers.Map();
        if (~isa(cfg.prop, 'containers.Map'))
            % use a non-empty synthetic key: Octave 8.4's containers.Map
            % crashes in sort_keys -> orderfields -> cell2struct when a key
            % is the empty string. The key never leaks: the single-key
            % short-circuit below unwraps the value before returning.
            prop('1') = cfg.prop;
        else
            prop = cfg.prop;
        end
        wavelengths = prop.keys;
        for waveid = prop.keys
            wv = waveid{1};
            pp = prop(wv);
            if (size(prop(wv), 1) < min([nn, ne])) % label based prop
                if (isfield(cfg, 'seg'))
                    if (length(cfg.seg) == nn)
                        bkprop(wv) = pp(cfg.seg(cfg.face(1)) + 1, :);
                    elseif (length(cfg.seg) == ne)
                        [xi, yi] = find(cfg.elem == cfg.face(1));
                        bkprop(wv) = pp(cfg.seg(xi(1)) + 1, :);
                    else
                        error('cfg.seg must match the length of node or elem');
                    end
                else
                    error('labeled proper is defined, but cfg.seg is not given');
                end
            elseif (size(pp, 1) == nn) % node based prop
                bkprop(wv) = pp(cfg.face(1), :);
            elseif (size(pp, 1) == ne) % elem based prop
                [xi, yi] = find(cfg.elem == cfg.face(1));
                bkprop(wv) = pp(xi(1), :);
            else
                error('row number of cfg.prop is invalid');
            end
        end
        if (length(bkprop.keys) == 1)
            bkprop = bkprop(wavelengths{1});
        end
    end
else
    if (isfield(cfg.bulk, 'mua'))
        bkprop(1) = cfg.bulk.mua;
    end
    if (isfield(cfg.bulk, 'dcoeff'))
        bkprop(2) = 1 / (3 * cfg.bulk.dcoeff);
        bkprop(3) = 0;
    end
    if (isfield(cfg.bulk, 'musp'))
        bkprop(2) = cfg.bulk.musp;
        bkprop(3) = 0;
    end
    if (isfield(cfg.bulk, 'g'))
        bkprop(3) = cfg.bulk.g;
    end
    if (isfield(cfg.bulk, 'n'))
        bkprop(4) = cfg.bulk.n;
    end
end
