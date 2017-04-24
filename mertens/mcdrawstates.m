function states = mcdrawstates(P, N, initndx)
% simulates evolution of Markov States for N steps starting from P(inindx,:)
% optional: initndx is index of inital state [default: unconditional distribution]
% function states = mcdrawstates(P, N, initndx)

error(nargchk(2,3,nargin))
if nargin < 3;
   p0 = mclimit(P);
elseif isscalar(initndx)
   p0 = P(initndx,:)';
else
   if initndx == 0
      p0 = mclimit(P);
   else
      p0 = initndx;
   end
end

states      = zeros(N, 1);
draws       = rand(N, 1);
states(1)   = find(draws(1) < cumsum(p0), 1, 'first');
Pcdf        = cumsum(P, 2);
for t = 2 : N
   states(t) = find(draws(t) < Pcdf(states(t-1), :), 1, 'first');
end