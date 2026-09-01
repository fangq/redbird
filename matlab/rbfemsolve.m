function varargout = rbfemsolve(Amat, rhs, method, varargin)
%
% res=rbfemsolve(Amat, rhs, method, options ...)
%
% Solving a linear system defined by Amat*res=rhs using various methods
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     Amat: the left-hand-size matrix, can be sparse
%     rhs:  the right-hand-size vector or matrix (multiple rhs)
%     method: (optional) a string specifying the solving method
%            'mldivide': use the left-divide method: res=Amat\rhs
%            'blqmr': use the block-QMR iterative method for entire RHS
%
%          other supported methods include
%            qmr, tfqmr, cgs, gmres, pcg, cgs, minres, symmlq, bicgstab;
%
%          if method is a positive integer, it calls blqmr.m (part of the
%          blit toolbox) with block size defined by "method".
%
%          for all non-block solvers, adding "par" prefix calls parfor to
%          solve each RHS in parallel, one must call matlabpool or parpool
%          (matlab 2016 or newer) to create the pool first
%     options: (optional) additional parameters are passed to the solver
%          function to specify solver options, please help qmr, cgs, etc to
%          see the additionally accepted parameters
% output:
%
%
% license:
%     GPL version 3, see LICENSE_GPLv3.txt files for details
%
% -- this function is part of Redbird toolbox
%

if (nargin < 3)
    method = 'mldivide';
end
if (nargout < 1)
    error('output can not be empty');
end
outputnum = nargout;

if (~ischar(method) && nargout > 0)
    res = zeros(size(rhs));
    if (method < 0)
        method = -method;
        len = 1:method:size(rhs, 2);
        res = cell(1, length(len) - 1);
        parfor i = 1:length(res)
            vo = cell(1, outputnum);
            maxidx = min([size(rhs, 2), len(i) + method - 1]);
            [vo{:}] = blqmr(Amat, full(rhs(:, len(i):maxidx)), varargin{:});
            res{i} = vo{1};
        end
        res = cell2mat(res);
    else
        vo = cell(1, nargout);
        for i = 1:method:size(rhs, 2)
            maxidx = min([size(rhs, 2), i + method - 1]);
            [vo{:}] = blqmr(Amat, full(rhs(:, i:maxidx)), varargin{:});
            res(:, i:maxidx) = vo{1};
        end
    end
    vo{1} = res;
    varargout = vo;
    return
end

% block solvers can handle multiple RHSs
switch method
    case 'blqmr'
        [varargout{1:nargout}] = blqmr(Amat, full(rhs), varargin{:});
        return
    case 'mldivide'
        varargout{1} = full(Amat \ rhs);
        for i = 2:nargout
            varargout{i} = [];
        end
        return
end

% non-block solvers have to solve one RHS at a time

sol = cell(1, size(rhs, 2));
for i = 1:size(rhs, 2)
    sol{i} = cell(1, nargout);
end

if (regexp(method, '^par'))
    method = regexprep(method, '^par', '');
    parfor i = 1:size(rhs, 2)
        switch method
            case 'qmr'
                [sol{i}{:}] = qmr(Amat, rhs(:, i), varargin{:});
            case 'tfqmr'
                [sol{i}{:}] = tfqmr(Amat, rhs(:, i), varargin{:});
            case 'cgs'
                [sol{i}{:}] = cgs(Amat, rhs(:, i), varargin{:});
            case 'lsqr'
                [sol{i}{:}] = lsqr(Amat, rhs(:, i), varargin{:});
            case 'gmres'
                [sol{i}{:}] = gmres(Amat, rhs(:, i), varargin{:});
            case 'gmres'
                [sol{i}{:}] = gmres(Amat, rhs(:, i), varargin{:});
            case 'pcg'
                [sol{i}{:}] = pcg(Amat, rhs(:, i), varargin{:});
            case 'minres'
                [sol{i}{:}] = minres(Amat, rhs(:, i), varargin{:});
            case 'symmlq'
                [sol{i}{:}] = symmlq(Amat, rhs(:, i), varargin{:});
            case 'bicgstab'
                [sol{i}{:}] = bicgstab(Amat, rhs(:, i), varargin{:});
            otherwise
                error(['method "', method, '" is not supported']);
        end
    end
else
    % do not use parfor
    for i = 1:size(rhs, 2)
        switch method
            case 'qmr'
                [sol{i}{:}] = qmr(Amat, rhs(:, i), varargin{:});
            case 'tfqmr'
                [sol{i}{:}] = tfqmr(Amat, rhs(:, i), varargin{:});
            case 'cgs'
                [sol{i}{:}] = cgs(Amat, rhs(:, i), varargin{:});
            case 'lsqr'
                [sol{i}{:}] = lsqr(Amat, rhs(:, i), varargin{:});
            case 'gmres'
                [sol{i}{:}] = gmres(Amat, rhs(:, i), varargin{:});
            case 'gmres'
                [sol{i}{:}] = gmres(Amat, rhs(:, i), varargin{:});
            case 'pcg'
                [sol{i}{:}] = pcg(Amat, rhs(:, i), varargin{:});
            case 'minres'
                [sol{i}{:}] = minres(Amat, rhs(:, i), varargin{:});
            case 'symmlq'
                [sol{i}{:}] = symmlq(Amat, rhs(:, i), varargin{:});
            case 'bicgstab'
                [sol{i}{:}] = bicgstab(Amat, rhs(:, i), varargin{:});
            otherwise
                error(['method "', method, '" is not supported']);
        end
    end
end

res = cell(1, nargout);
for i = 1:nargout
    try
        res{i} = cell2mat(cellfun(@(x) x{i}, sol, 'uniformoutput', true));
    catch
        res{i} = cell2mat(cellfun(@(x) x{i}, sol, 'uniformoutput', false));
    end
end

varargout = res;
