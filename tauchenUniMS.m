%% MSVAR Tauchen (1986): Uniform grid
%Requires: Neurel Network Toolbox -- combvec.m
clear all;
%% Discretization Paramters
nS = 2; %nS
n = 2; %Number of elements in VAR
nPoints = NaN(n,nS);
nPoints = [5,5; 5,5]; %Number of points of variable j in regime i
nSig = 3.25; %Number of std. dev. from mean to make grid;

%% MS-VAR(1) Specification
Pi = [0.97,0.03; 0.95,0.05];
PhiCell = cell(nS,1);
muCell = cell(nS,1);
SigCell = cell(nS,1);
muCell = cell(nS,1);
cCell = cell(nS,1);
PhiCell{1} = [0.7,0.1; 0,0.9];
PhiCell{2} = [0.1,0;0,0.99];
muCell{1} = [0.002;0.003];
muCell{2} = zeros(n,1);
SigCell{1} = sqrtm([0.1,0; 0,0.3]);
SigCell{2} = sqrtm([0.1,0; 0,0.1]);
vCell{1} = [0.1,0; 0,0.3];
vCell{2} = [0.1,0; 0,0.1];
cCell{1} = (eye(n)-PhiCell{1})\muCell{1};
cCell{2} = (eye(n)-PhiCell{2})\muCell{2};

%Check for MSS - following Cho (2016, RED) - Doesn't make senSe to discretize
    %unbounded process
    Rss = [Pi(1,1).*kron(PhiCell{1},PhiCell{1}),...
        Pi(2,1).*kron(PhiCell{1},PhiCell{1});...
        Pi(1,2).*kron(PhiCell{2},PhiCell{2}),...
        Pi(2,2).*kron(PhiCell{2},PhiCell{2})];

    if max(eig(Rss))<1
        disp('MS-VAR is MSS')
    else
        error('MS-VAR specification not MSS as defined in Cho (2016)')
    end

%Get ergodic distribution of Markov chain
q  = getStatMarkov(Pi);

%Calculate Unconditional Mean and Variance
tSum = zeros(n,1);
for ii=1:nS
    tSum = tSum + q(ii).*muCell{ii};
end
Ex = tSum;

varCell = cell(nS,1);
for ii=1:nS
    varCell{ii} = ...
        reshape((eye(n^2)-kron(PhiCell{ii},PhiCell{ii}))\vCell{ii}(:),n,n);
end

tSum = zeros(n);
for ii=1:nS
    tSum = tSum + q(ii).*varCell{ii};
end
VarX  = tSum;


%Create \Pi_{i->j}
Pi_ij = cell(nS);
for ii=1:nS
    for jj=1:nS
       Pi_ij{ii,jj} = fn_var_to_markov(eye(n),muCell{jj},PhiCell{jj},...
           vCell{jj},nPoints(:,jj),1000,1);
    end
end

%Create \Pi from \Pi_{i->j}
tSize = 0;
sizeMat = NaN(nS,1);
for ii=1:nS
    tSize = tSize + prod(nPoints(:,ii));
    sizeMat(ii) =  prod(nPoints(:,ii));
end
endIII = cumsum(sizeMat);
nMC = tSize;
initIII = NaN(nS,1);
initIII(1)=1;
initIII(2:end) = endIII(1:end-1)+1;


outPi = NaN(nMC);
cnt = 1;
for ii=1:nS
    for iii=1:sizeMat(ii)
        for jj=1:nS
            outPi(cnt,initIII(jj):endIII(jj)) = Pi(ii,jj)*...
                Pi_ij{ii,jj}(iii,:);
        end
        cnt = cnt + 1;
    end
end
