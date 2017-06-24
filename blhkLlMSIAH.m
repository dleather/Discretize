function [ll,fp,sp] = blhkLlMSIAH(data,Sig,A,v,P,smooth)
%{
This code follows Krolzig (1997) in writing code to efficiently compile the
log-likelihood of the the MS-VAR(1) model with switching intercepts,
heteroskedasticity, and switching autoregressive coefficients.

IN
    data :: (T* x K) time series data where T is length of TS and K is #
        variables in VAR
    Sig :: (K x KM) matrix of covariance matrices stackes horizontally. M
        is number of total regimes (compound if multiple Markov Chains)
    A :: (K x KM) matrix of autoregressive matrices
    v :: (K x M) matrix of intercepts
    P :: (M x M) transition matrix where rows sum to unity
    smooth :: dummy to produce smooth probabiltiies (backwards)

OUT
    ll :: value of log-likelihood
    fp :: (T x M) matrix of filtered probabilities
    sp :: (T x M) matrix of smoothed probabilities (following Kim (1994)   
        method)
%}

%% Input Check
if nargin <6
    smooth=0;
    sp = [];
end
[Tdata,K] = size(data);
[M,tM] = size(P); if M~=tM, error('P is not square'); end
[tK,KM] = size(Sig); 
if (K~=(KM/M))||(tK~=K), error('Sig/P or Sig/data dimension mismatch');end
if max(size(A)~=[K,KM]), error('A/Sig dimension mismatch'); end
if max(size(v)~=[K,M]), error('v dimension mismatch'); end
if ~((smooth==0)||(smooth==1)), error('smooth is not dummy'); end
if smooth==0, sp = NaN; end


%% Set-Up State Space - See Table 9.11 Krozlig (1997)
% Measurement Equation: y_t = (\xi'_t \otimes X_t)vec(B) + u_t
% State or xition equation: \xi_{t+1} = F \xi_t + v_{t+1}

T = Tdata-1;
Y = data(2:end,:);
F = P';
barX = [ones(T,1),data(1:end-1,:)];
y = Y'; y = y(:);

tSig = Sig;
Sig = kron(eye(M),ones(K));
Sig(logical(Sig)) = tSig(:);

B = NaN(K*(1+K),M);

for m=1:M
   B(:,m) = [v(:,m); reshape(A(1:K,K*(m-1)+1:K*m),K*M,1)];

end    


%% Form Likelihood
%Get ergodic distribution
q = (approxStatMarkov(P,1e-06,100000))';

%preallocate
xiHat_ttm1 = NaN(M,T+1);
AA = NaN(K,M*K);
srp2fac = (2*pi)^(-K/2);
eta = NaN(M,T);
eyeM = eye(M);
oneM = ones(1,M);
lhMat = NaN(T,1);


%Get (T x KM) matrix of innovations, e
YY = repmat(Y,1,M);
vv = reshape(v,1,K*M);
for m=1:M
    AA(1:K,K*(m-1)+1:K*m) = A(1:K,K*(m-1)+1:K*m)';
end
BB = [vv;AA];
e = YY - barX*BB;

%Initiate recursion
xiHat_ttm1(:,1) = q;


%Form likelihood
for t=1:T
    %get eta
    for m=1:M
       sig = tSig*(kron(eyeM(:,m),eye(K))); %Pick out Sigma
       te = e(t,K*(m-1)+1:K*m);
       eta(m,t) = srp2fac.*(det(sig)^(-0.5)).*exp(-0.5*(te/sig*te'));
    end
    
    %update filtered probability
    tMat = (eta(:,t).*xiHat_ttm1(:,t));
    xiHat_ttm1(:,t+1) = (F*tMat)/(oneM*tMat);
    
    %store likelihood
    lhMat(t) = eta(:,t)'*xiHat_ttm1(:,t); 
end

fp = xiHat_ttm1(:,2:end)';
ll = sum(log(lhMat));

if smooth==1
   %preallocate
   xiHat_tT = NaN(M,T);
   
   %initiate recursion
   xiHat_tT(:,T) = F*xiHat_ttm1(:,end-1);
   
   %loop
   for t=T-1:-1:1
       xiHat_tT(:,t) = (F'*(xiHat_tT(:,t+1)./xiHat_ttm1(:,t+1))).*...
           (F*xiHat_ttm1(:,t));
       xiHat_tT(:,t) = xiHat_tT(:,t)./sum(xiHat_tT(:,t));
       
   end
   
    sp = xiHat_tT';
end

%% Output Check
if ~(isscalar(ll)), error('ll output is not scalar'); end
if any((size(fp)==[T,M])==0), error('fp dimensions wrong'); end
if smooth==1
    if any((size(sp)==[T,M])==0), error('fp dimensions wrong'); end
end

end




