function [ Pi ] = r2pmatCmb(theta)

%{
r2pmatCmb contructs a (2^(N/2) x 2^(N/2)) Markov transition matrix from (N x 1)
vector theta \in \mathbb{R}^(N/2). The idea is to use the normalization used in
Ding (2012) in order to map  f: R --> [0,1] s.t. f(x) = 1 / (1 + e^x)
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

tTheta = reshape(theta,2,nTheta/2);
tFn = @(x) 1 / ( 1 + exp(x));
%Construct 2x2 matrices
piCell = cell(nTheta/2,1);
for ii=1:(nTheta/2)
   piCell{ii} = [tFn(tTheta(1,ii)),1-tFn(tTheta(1,ii));...
       1-tFn(tTheta(2,ii)), tFn(tTheta(2,ii))];
end

if (nTheta/2)>1
%Combine 2x2 matrices 
temp = kron(piCell{1},piCell{2});
if (nTheta/2)>2
    for ii=3:(nTheta/2)
        temp = kron(temp,piCell{3});
    end
end

Pi = temp;
else
    Pi = piCell{1};
end

%Output check
if ~isequal(sum(Pi,2),ones(2^(nTheta/2),1))
    error('Rows don''t sum to unity')
end

end
