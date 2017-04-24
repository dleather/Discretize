function [flag, msg, absorbing] = istransitionmat(P, tol)
% checks whether matrix is a valid transitionmatrix

error(nargchk(1,2,nargin))

if nargin < 2
   tol = 1e-12;
end

msg   = [];
flag  = true;

[Ns, tmp] = size(P);
if Ns ~= tmp
   flag = false;
   msg  = sprintf('%s not square.', msg);
end

if ~all(abs(sum(P, 2) - 1) < tol)
   flag = false;
   msg  = sprintf('%s rows do not sum to unity.', msg);
end

if nargout > 2
   absorbing = any(P(:) == 0);
end