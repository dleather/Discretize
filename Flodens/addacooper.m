function [Z,PI] = addacooper(n,mu,rho,sigma);
% Approximate n-state AR(1) process following Tauchen (1986) and Tauchen & Hussey (1991). 
% See Adda & Cooper (2003) pp 57-.
%
% Z(t+1) = mu*(1-rho) + rho*Z(t) + eps(t+1)
%
% where  std(eps) = sigma
%
% Martin Flodén, 2005

sigmaUNC = sigma/sqrt(1-rho^2);
E  = zeros(n+1,1);
Z  = zeros(n,1);
PI = zeros(n,n);
MFPI = zeros(n,n);

E(1)   = -1E6;
E(end) = 1E6;
for i = 2:n
    E(i) = sigmaUNC*norm_inv((i-1)/n) + mu;
end

for i = 1:n
    Z(i) = n*sigmaUNC*(norm_pdf((E(i)-mu)/sigmaUNC) - norm_pdf((E(i+1)-mu)/sigmaUNC)) + mu;       
end

for i = 1:n
    for j = 1:n
        E1 = E(j);
        E2 = E(j+1);
        th_fcn  = @(u) n/sqrt(2*pi*sigmaUNC^2) * (exp(-(u'-mu).^2 / (2*sigmaUNC^2)) .* ...
                      (norm_cdf((E2-mu*(1-rho)-rho*u')/sigma) - norm_cdf((E1-mu*(1-rho)-rho*u')/sigma)));
        PI(i,j) = quadl(th_fcn,E(i),E(i+1),1e-10);
        MFPI(i,j) = norm_cdf((E(j+1)-mu*(1-rho)-rho*Z(i))/sigma) - norm_cdf((E(j)-mu*(1-rho)-rho*Z(i))/sigma);
       
    end
end

for i = 1:n
    PI(i,:) = PI(i,:) / sum(PI(i,:));
    MFPI(i,:) = MFPI(i,:) / sum(MFPI(i,:));
end


function c = norm_cdf(x)
    c = 0.5 * erfc(-x/sqrt(2));

function p = norm_pdf(x)
    p = 1/sqrt(2*pi) * exp(-(x)^2/2);

function y = norm_inv(x)
    f = @(x0,a) (0.5 * erfc(-x0/sqrt(2))-a);
    y = fzero(f,0,[],x);
        