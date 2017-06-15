function [ Pi ] = r2pmat(theta)

%{
r2pmat follows the algorithm described in Ding (2012)
"An Implementation of Markov Regime Switching Model with Time Varying Transition
 Probabilities in Matlab", found here:
 https://www.researchgate.net/publication/256021580_An_Implementation_of_Markov_Regime_Switching_Model_with_Time_Varying_Transition_Probabilities_in_Matlab
 and as implemented without state dependence in Perlin's MSRegress package
 found https://github.com/msperlin/MS_Regress-Matlab .

 The goal is to construct a (kxk) Markov transition matrix whose rows sum
 to unity, from k(k-1) free parameters on the real number line. This transforms
 the problem of estimation a MS-VAR from a constrained optimization problem to
 an unconstrained optimization problem.

 IN

    theta :: (k(k-1) x 1) vector in R

 OUT

    Pi :: (k x k) Markov transition matrix whose rows sum to unity

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

%Get k
k  = max((1-sqrt(1+4*nTheta))/2,(1+sqrt(1+4*nTheta))/2);
if k~=round(k)
    error('N is not an integer')
end

theta = reshape(theta,k-1,k);

p=ones(k);
temp_p=ones(k);
for i=1:k-1
    for j=1:k
        temp_p(i,j)=normcdf(theta(i,j));
        p(i+1,j)=p(i,j)*(1-temp_p(i,j));
    end
end

Pi = (p.*temp_p)';

%Output check
if ~isequal(sum(Pi,2),ones(k,1))
    error('Rows don''t sum to unity')
end

end
