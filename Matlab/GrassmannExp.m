function U1 = GrassmannExp(U0,Delta)
% Inputs:
%   U0: Stiefel representative of subspace
%   Delta: Tangent vector in horizontal space at U0
% Output:
%   Uend: Grassmann exponential of Delta at U0

[Q,Sigma,V] = svd(Delta,0);
cosSigma = diag(cos(diag(Sigma)));
sinSigma = diag(sin(diag(Sigma)));

U1 = (U0*V*cosSigma + Q*sinSigma)*V';

return
end