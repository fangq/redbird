function [finalAmat, finalrhs, nblock] = rbcreateinv(Amat, ymeas, ymodel, params, output)
%
% [finalAmat, finalrhs, nblock]=rbcreateinv(Amat, ymeas, ymodel, params, output)
%
% Convert an inverse problem to the log-amplitude and unwrapped phase form
% i.e. Amat*x=[ymeas - ymodel] =>
% [Alogamp Aphase]*x=[log10(ymeas)-log10(ymodel);angle(ymeas)-angle(ymodel)]
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     Amat: a complex-valued LHS matrix
%     ymeas: complex valued measurement data
%     ymodel: complex valued model prediction
%     param: unknown parameters for multiple wavelengths
%     output: 'complex': complex-valued system,
%                     finalAmat=Amat; finalrhs=ymeas-ymodel
%             'real': real-valued system,
%                     finalAmat=[Re(Amat) -Im(Amat);Im(Amat) Re(Amat)]
%                     finalrhs =[Re(ymeas-ymodel) ; Im(ymeas-ymodel)]
%             'logphase': real-valued system,
%                     finalAmat=see manual
%                     finalrhs =[log10(ymeas)-log10(ymodel) ; angle(ymeas)-angle(ymodel)]
%
% output:
%     finalAmat: converted LHS matrix, see above
%     finalrhs: converted RHS matrix, see above
%     nblock: after conversion, x can be converted to [Re(x);Im(x)]; if
%          this happens, nblock is set to 2, otherwise 1
%
% license:
%     GPL version 3, see LICENSE_GPLv3.txt files for details
%
% -- this function is part of Redbird-m toolbox
%

% convert an inverse problem to the log-amplitude and unwrapped phase form
%

if (nargin < 5)
    output = 'complex';
end

if (nargin < 4)
    param = [];
end

wavelengths = {'_'};

if (isa(Amat, 'containers.Map'))
    wavelengths = Amat.keys;
else
    Amat = containers.Map({'_'}, {Amat});
    ymeas = containers.Map({'_'}, {ymeas});
    ymodel = containers.Map({'_'}, {ymodel});
end

newA = containers.Map();
newrhs = containers.Map();

nblock = 1;

if (strcmp(output, 'complex'))
    newA = Amat;
    for waveid = wavelengths
        wv = waveid{1};
        newrhs(wv) = ymeas(wv) - ymodel(wv);
    end
else
    for waveid = wavelengths
        wv = waveid{1};
        rhs = ymeas(wv) - ymodel(wv);
        if (strcmp(output, 'real'))
            newA(wv) = real(Amat(wv));
            newrhs(wv) = real(rhs);

            if (~isreal(rhs) && ~isreal(Amat(wv)))
                newA(wv) = [real(Amat(wv)) -imag(Amat(wv))
                            imag(Amat(wv)) real(Amat(wv))];
                newrhs(wv) = [real(rhs); imag(rhs)];
                nblock = 2;
            end
        elseif (strcmp(output, 'logphase'))
            temp = conj(ymodel(wv)) ./ abs(ymodel(wv) .* ymodel(wv));
            temp = repmat(temp(:), 1, size(Amat(wv), 2)) .* Amat(wv);
            if (isreal(ymodel(wv)))
                newA(wv) = real(temp);
                temp = log(abs(ymeas(wv))) - log(abs(ymodel(wv)));
                newrhs(wv) = temp(:);
            else
                newA(wv) = [real(temp); imag(temp)];
                temp = log(abs(ymeas(wv))) - log(abs(ymodel(wv)));
                ptemp = angle(ymeas(wv))   - angle(ymodel(wv));
                newrhs(wv) = [temp(:); ptemp(:)];
                nblock = 2;
            end
        end
    end
end

finalAmat = [];
finalrhs = [];

if (isa(params, 'struct'))
    paramlist = intersect(fieldnames(params), {'hbo', 'hbr', 'water', 'lipids', 'aa3'});
    if (isempty(paramlist))
        error('specified parameters are not supported');
    end
    extins = rbextinction(wavelengths, paramlist);
    for i = 1:length(wavelengths)
        wv = waveid{i};
        wvAmat = [];
        for j = 1:length(paramlist)
            wvAmat = [wvAmat, newA(wv) * extins(i, j)];
        end
        finalAmat = [finalAmat; wvAmat];
        finalrhs = [finalrhs; newrhs(wv)];
    end
else
    for i = 1:length(wavelengths)
        wv = waveid{i};
        finalAmat = [finalAmat; newA(wv)];
        finalrhs = [finalrhs; newrhs(wv)];
    end
end
