function res = rbreginvover(Amat, rhs, lambda, LTL, blocks, varargin)
%
% res=rbreginvover(Amat, rhs, lambda, LTL, blocks, varargin)
%
% Solve an overdetermined Gauss-Newton normal equation
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     Amat: the left-hand-side matrices (a containers.Map object) at specified wavelengths
%     rhs: the right-hand-side vectors for all sources (independent of wavelengths)
%     lambda: the Tikhonov regularization parameter
%     LTL: the regularization matrix in the form of the inverse of the
%         unknown-covariance (inv(C_mu)), which is L'*L
%     blocks: the dimensions of the 2D submatrices of Amat, needed if
%         size(LTL,1)==size(LTL,2) is not the same as size(Amat,2)
%
% output:
%     res: the least-square solution of the matrix equation
%
% license:
%     GPL version 3, see LICENSE_GPLv3.txt files for details
%
% -- this function is part of Redbird toolbox
%

% solve an overdetermined Gauss-Newton normal equation
%  delta_mu=inv(J'J + lambda*(L'L))*J'*(y-phi)

% if any node has no sensitivity, remove them from inversion
length0 = size(Amat, 2);
idx0 = find(sum(Amat) ~= 0);
if (length(idx0) < length0)
    Amat = Amat(:, idx0);
    % TODO: need to shrink Lqr as well
end

emptydata = find(sum(Amat') ~= 0);
if (length(emptydata) < size(Amat, 1))
    Amat = Amat(emptydata, :);
    rhs = rhs(emptydata);
end

if (nargin >= 4 && ~isempty(LTL))
    nx = size(LTL, 1);
    oldnx = length0 / length(fieldnames(blocks));
    if (nx > length(idx0) / length(fieldnames(blocks)))
        Lidx = idx0(idx0 <= nx);
        LTL = LTL(Lidx, Lidx);
        oldnx = nx;
        nx = size(LTL, 1);
    end
end

rhs = Amat' * rhs(:);

Hess = Amat' * Amat; % Gauss-Hessian matrix, approximation to Hessian (2nd order)

% [Hess,Gdiag]=rbnormalizediag(Hess);

if (nargin < 4 || isempty(LTL))
    Hess(1:1 + size(Hess, 1):end) = Hess(1:1 + size(Hess, 1):end) + lambda;
else
    if (size(Hess, 1) == size(LTL, 1))
        Hess = Hess + lambda * LTL;
    else
        nx = size(LTL, 1);
        for i = 1:nx:size(Hess, 1)
            Hess(i:i + nx - 1, i:i + nx - 1) = Hess(i:i + nx - 1, i:i + nx - 1) + lambda * LTL;
        end
        %         if(isempty(blocks)) % assume the Hess matrix size is multiples of LTL
        %             for i=1:nx:size(Hess,1)
        %                 Hess(i:i+nx-1,i:i+nx-1)=Hess(i:i+nx-1,i:i+nx-1)+lambda*LTL;
        %             end
        %         else
        %             len=cumsum([1; structfun(@(x) x(2), blocks(1))]);
        %             for i=1:length(blocks)
        %                 if(nx==len(i+1)-len(i)+1)  % if the block size match LTL
        %                     Hess(len(i):len(i+1)-1,len(i):len(i+1)-1)=Hess(len(i):len(i+1)-1,len(i):len(i+1)-1)+lambda*LTL;
        %                 else % if size does not match, add lambda*I
        %                     idx=sub2ind(size(Hess),len(i):len(i+1)-1,len(i):len(i+1)-1);
        %                     Hess(idx)=Hess(idx)+lambda;
        %                 end
        %             end
        %         end
    end
end
[Hess, Gdiag] = rbnormalizediag(Hess);
res = Gdiag(:) .* rbfemsolve(Hess, Gdiag(:) .* rhs, varargin{:});
% res = rbfemsolve(Hess,rhs,varargin{:});

if (length(idx0) < length0)
    res0 = zeros(length0, size(res, 2));
    res0(idx0, :) = res;
    res = res0;
end
