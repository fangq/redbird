function lambda = rbregemperical(Hess, residual, alpha)
%
% lambda=rbregemperical(Hess, residual, alpha)
%
% Apply an emperical approach to determine the regularization coeff lambda
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     Hess: the Gauss-Hessian matrix
%     residual: the total data-model misfit (residual) of the previous iteration
%     alpha: an emperical regularization parameter
%
% output:
%     lambda: the estimated emperical regularization parameter
%
% license:
%     GPL version 3, see LICENSE_GPLv3.txt files for details
%
% -- this function is part of Redbird toolbox
%

ggav = mean(diag(Hess));
lambda = alpha * ggav * residual * residual;
