%PINVF: Pseudoinverse based on r singular values of A.
%	X = PINVF(A,r,tol) produces a matrix X of the same dimensions
%	as A' so that A*X*A = A, X*A*X = X and AX and XA
%	are Hermitian. The computation is based on SVD(A) and r
%	singular values.
%
%       X = pinvf(A,r,tol)
%
%	See also RANK.
%
%       L.G. Van Willigenburg, W.L. De Koning, 28-11-95.
%
  function X = pinvf(A,r,tol)

  [U,S,V] = svd(A,0);
  if min(size(S)) == 1
     S = S(1);
  else
     S = diag(S);
  end
  if nargin~=3
    tol = max(size(A)) * S(1) * eps;
  end
  rt = sum(S > tol);
  if nargin>1
     rt=min(rt,r);
  end
  if (rt == 0)
     X = zeros(size(A'));
  else
     S = diag(ones(rt,1)./S(1:rt));
     X = V(:,1:rt)*S*U(:,1:rt)';
  end

