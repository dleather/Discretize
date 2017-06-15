function [cov,sig] = sphere2cov(theta)
%{
Following Pinhiero & Bates, "Unconstrained Parameterizations for
Variance-Covariance Matrices", this function takes N(N+1)/2 inputs, \theta_i,
where N is the dimension of the square var-cov matrix. It outputs both the Cov
matrix and cholesky decomposition following the "Spherical Parameterization"

cov = sig' * sig

IN
    theta :: (N(N+1)/2 x 1) vector of unconstrained parameters. First
             N correspond to the diagonals of sig, and latter N(N-1)/2 refer to
             off-diagonals
OUT
    cov = sig' * sig :: (N x N) covariance matrix ::
    sig :: cholesky decomposition of cov

%}

[nTheta,m] = size(theta);
if m>1
    if nTheta == 1
        theta = theta';
        nTheta = m;
    else
        error('theta is not vector')
    end
end

%Get N
N  = max((-1-sqrt(1+8*nTheta))/2,(-1+sqrt(1+8*nTheta))/2);
if N~=round(N)
    error('N is not an integer')
end

%Get matrix l descrbed below
l = zeros(N); %Matrix where [l_i]_j corresponds to j,i element of l
for ii=1:N
    l(ii,1) = exp(theta(ii));
end

for ii=2:N
    for jj=2:ii
        l(ii,jj) = (exp(theta(N+(ii-2)*(ii-1)/2+(jj-1)))*pi)/...
            (1+exp(theta(N+(ii-2)*(ii-1)/2+(jj-1))));
    end
end

%Construct L
L = zeros(N);
for ii=1:N
    if ii>1
        for jj=ii:-1:1
            if (jj<ii)&&(jj>1)
                L(jj,ii) = l(ii,1) * prod(sin(l(ii,2:jj)))*cos(l(ii,jj+1));
            elseif jj==ii
                L(jj,ii) = l(ii,1) * prod(sin(l(ii,2:jj)));
            else %jj==1
                L(jj,ii) = l(ii,1) * cos(l(ii,2));
            end
        end
    else
        L(ii,ii) = l(ii,ii);
    end
end

sig = L';
cov = L'*L;

end
