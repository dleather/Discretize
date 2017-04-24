function [p0, iterations] = mclimit(P, iterflag, tol, maxiter)
% limting dstribution of markovchain
% function [p0, iterations] = mclimit(P, iterflag, tol, maxiter)
% iterflag = true: pick iterative solution
% iterflag = false: look for unit eigenvector
% leave iterflag empty to let mclimit pick the efficient algorithm (iterative for Ns > 100)


error(nargchk(1,3,nargin))


if nargin < 4
   maxiter = 1000;
end

if nargin < 3
   tol = 1e-12;
end

if nargin < 2
   iterflag = [];
end


if any(P(:) == 0)
   warning('elmi:base', 'absorbing states!')
end

Ns = length(P);

if isempty(iterflag)
   if Ns > 100
      iterflag = true;
   else
      iterflag = false;
   end
end

if Ns > 100 && ~iterflag
   warning('elmi:base', 'Eigenvalue computation slow with %d states', Ns);
end

if iterflag
   Pold  = P^20;
   delta = tol * 100;
   iterations  = 1;
   while delta > tol && iterations < maxiter
      Pnew        = Pold^5;
      delta       = sum(abs(Pnew(:) - Pold(:)));
      iterations  = iterations + 1;
      Pold        = Pnew;
   end

   if iterations >= maxiter
      warning('elmi:base', 'Termination because Maxiter %d reached', maxiter);
   end
   
   p0 = Pnew(1,:)';
   
   % check that rows are really equal
   if mae(p0', mean(Pnew)) > tol || mean(std(Pnew)) > tol
      error('there are still differences between the rows')
   end
   
else % seek unit eigenvalue
   [V, D] = eig(P');
   lambdas = diag(D);
   ndx = lambdas == 1;
   if all(ndx == 0)
      ndx = abs(lambdas - 1) < tol;
      if all(ndx == 0)
         error('Cannot find a unit root in transition matrix')
      end
   end
   p0 = V(:, ndx);
   p0 = p0 / sum(p0);
end