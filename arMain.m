%{
SCRIPT: Approximate AR(1) 
%}

%% Model Parameters: y_t = \phi y_t-1 + \sig \varepsilon_t
phi = 0.9; %AR coefficient
sig = 0.01; %std. dev of shock
barX = 0;
N = 20; %Number of nodes
%omega = @(x) normpdf(x); % Set omega to std. normal PDF


%% Construct f(\bar(y)_k | x)

%get nodes and weights via gauss-quadrature
[u,w] = hernodes(N);
xNodes = sqrt(2)*sig*u + barX;
piMat = NaN(N); %Markov matrix for discretized process
for ii=1:N
    tMat = NaN(N,1);
    for jj=1:N
        tMat(jj) = exp(-((u(jj)-phi*u(ii)-barX)^2 -(u(jj)-barX)^2))...
            *w(jj)/sqrt(pi);
    end
    sX = sum(tMat);
    piMat(ii,:) = tMat'/sX;
end

%% Check via simulation
nSim = 1000;
nTs = 1000;

tsCell = cell(nSim);
%Simulate nSim series of nTs
for ii=1:nSim
    tsCell{ii} = simMC(xNodes,piMat,0,nTs,N);
end

bMat = NaN(nSim,1);
for ii=1:nSim
    YY = tsCell{ii}(2:end);
    XX = tsCell{ii}(1:end-1);
    bMat(ii) = (XX'*XX)\(XX'*YY);
end

[mean(bMat),prctile(bMat,5),prctile(bMat,95)]

%% simulate on continuous data
tsCCell = cell(nSim);
mdl = arima('ARLags',1,'Constant',0,'Variance',sig^2,'AR',phi);
for ii=1:nSim
   tsCCell{ii} = simulate(mdl,nTs); 
end

bbMat = NaN(nSim,1);
for ii=1:nSim
    YY = tsCCell{ii}(2:end);
    XX = tsCCell{ii}(1:end-1);
    bbMat(ii) = (XX'*XX)\(XX'*YY);
end
[mean(bbMat),prctile(bbMat,5),prctile(bbMat,95)]


function [ts] = simMC(u,Pi,initInd,nTs,N)
    %u :: Nx1 set of nodes
    %Pi :: NXN Markov transition matrix
    %initInd :: {0,1,...,N} initial node
    
    S = NaN(nTs,1);
    ts = NaN(nTs,1);
    if initInd==0
        q = getStatMarkov(Pi);
        tmp = rand;
        tNdx = 1;
        trig=0;
        while (tNdx<=N)&&(trig==0)
           if tmp<=sum(q(1:tNdx))
               S(1) = tNdx;
               ts(1) = u(tNdx);
               trig = 1;
           else
               tNdx = tNdx + 1;
           end
        end
    end
    
    for tt=2:nTs
        tmp = rand;
        tNdx = 1;
        trig=0;
        while (tNdx<=N)&&(trig==0)
           if tmp<=sum(Pi(S(tt-1),1:tNdx))
               S(tt) = tNdx;
               ts(tt) = u(tNdx);
               trig = 1;
           else
               tNdx = tNdx + 1;
           end
        end
    end            
end


        
    
    
    


