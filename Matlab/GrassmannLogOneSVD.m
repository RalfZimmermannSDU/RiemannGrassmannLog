function [Delta,U1star] = GrassmannLogOneSVD(U0,U1)
% Input
%   U0,U1: Stiefel representatives of subspaces
% Output
%   U1star: Adapted Stiefel representative of U1
%   Delta: Tangent vector in horizontal space at U0, from U0 to U1star
%
% This version of the Grassmann logarithm avoids the computation of
% the SVD of "H = U1star-U0*(U0'*U1star)"
% A short consideration shows that all the required information are 
% already encoded in the SVD of M = U1'*U0
%
% This modification of the Grassmann logarithm of Alg. 5.3 
% corresponds to the remark in Section 5.2 of the paper.
%

% Step 1: Procrustes
M = U1'*U0;
[Q1,S1,R1] = svd(M);
% Reorder the columns
Q1 = flip(Q1,2);
R1 = flip(R1,2);
S1 = diag(flip(diag(S1)));

% Calculate new rep
U1star = U1*Q1;

% Step 3: SVD without actual SVD
H = U1star-U0*(U0'*U1star);
singvals  = sqrt(1-diag(S1).^2);

Q2 = H./singvals';

Sigma = diag(asin(singvals));

% Step 3: Tangent vector
Delta = Q2*Sigma*R1';

return
end
