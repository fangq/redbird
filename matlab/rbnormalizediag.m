function [A, di] = rbnormalizediag(A0)
%
% [A, di]=rbnormalizediag(A0)
%
% Normalize/equalize a square matrix A0 and make the diagnonal all 1s, so that
% A0=di'*A*di
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     A0: the square matrix to be normalized
%
% output:
%     A: the square mtarix after normalization with all diagonal elements 1
%
% license:
%     GPL version 3, see LICENSE_GPLv3.txt files for details
%
% -- this function is part of Redbird toolbox
%

Adiag = diag(A0);
di = 1 ./ sqrt(Adiag);
A = (di * di') .* A0;
