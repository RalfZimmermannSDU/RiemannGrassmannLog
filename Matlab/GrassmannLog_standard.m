function [Delta] = GrassmannLog_standard(U0,U1)
% Input
%   U0,U1: Stiefel representatives of subspaces
% Output
%   Delta: Tangent vector in horizontal space at U0, 
%          from U0 to subspace spanned by U1
%
% This implementation of the algorithm follows
%
% P.-A. Absil, R. Mahony, and R. Sepulchre. 
% "Riemannian geometry of Grassmann manifolds with 
%  a view on algorithmic computation."
% Acta Applicandae Mathematica, 80(2):199–220, 2004. 
% doi:%10.1023/B:ACAP.0000013855.14971.91.
%
%
%




% Step 1: (I-U0*U0')*U1*(U0'*U1)^-1
M = U0'*U1;
N = U1/M - U0;

% Step 2: SVD
[Q2,S2,R2] = svd(N,0);
Sigma = diag(atan(diag(S2)));

% Step 3: Tangent vector
Delta = Q2*Sigma*R2';

return
end