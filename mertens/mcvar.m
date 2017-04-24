function Omega = mcvar(P, y, init, k)
% computes ACF of markov chain for y with k lags
% function Omega = mcvar(P, y, init, k)
% note: if k = 0 and init a vector, P=[] is admissible

error(nargchk(2,4,nargin))

if nargin < 4
   k = 0;
end
if nargin < 3 
   init = 0;
end

% get initial probabilities
if isscalar(init)
   if init == 0
      p    = mclimit(P)';
   else
      p = P(init, :);
   end
else
   p = init(:)';
end

% deviations from mean
[Ns, Ny] = size(y);
e  = y - repmat(p * y, Ns, 1);
switch k
   case 0
      enext = e;
   case 1
      enext = P * e;
   otherwise
      enext = P^k * e;
end

% Omega
Omega = zeros(Ny);
eprob    = e .* repmat(p', 1, Ny); % this is for speed
for s = 1 : Ns
   Omega  = Omega + eprob(s,:)' * enext(s,:);
   % old version:
%    Omega  = Omega + p(s) * e(s,:)' * enext(s,:);
%    if mae(p(s) * e(s,:)' * enext(s,:), eprob(s,:)' * enext(s,:)) > eps
%       warning('check it out: %12.6f', mae(p(s) * e(s,:)' * enext(s,:), eprob(s,:)' * enext(s,:)))
%       keyboard
%    end
end