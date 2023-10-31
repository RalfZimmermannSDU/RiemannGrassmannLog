function [Delta,U1star] = GrassmannLog(U0,U1)
% Input
%   U0,U1: Stiefel representatives of subspaces
% Output
%   U1star: Adapted Stiefel representative of U1
%   Delta: Tangent vector in horizontal space at U0, from U0 to U1star
%
% This version of the Grassmann logarithm corresponds
% to Alg. 5.3 of the associated paper.
%
%
%

% Step 1: Procrustes
M = U1'*U0;
[Q1,~,R1] = svd(M);
U1star = U1*(Q1*R1');

% Step 2: SVD
H = U1star-U0*(U0'*U1star);
[Q2,S2,R2] = svd(H,0);
Sigma = diag(asin(diag(S2)));

% Step 3: Tangent vector
Delta = Q2*Sigma*R2';

return
end
