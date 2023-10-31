function [ dist ] = subspaceDist(U0, U1)
%-------------------------------------------------------------
% compute Grassmann distance dist(U0,U1)
% in terms of the norm of the canonical angles 
%-------------------------------------------------------------

[Q, S, R] = svd((U0'*U1));

% optional: catch numerical round off errors:
% enforce that all singular vals are <=1
s = diag(S);
% for k = 1:length(s)
%     %if abs(s(k) - 1.0) < 100*eps || s(k)>1
%     if s(k)>1 
%         s(k) =1.0;
%     end
% end
theta = real(acos(s));
dist = norm(diag(theta),'fro');
return;
end

