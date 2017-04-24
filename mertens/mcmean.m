function Ybar = mcmean(P, y)
% computes mean of markov chain for y
% function Ybar = mcmean(P, y)

if isvector(P)
   P = P(:)';
end
Ybar = P * y;