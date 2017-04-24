%{
SCRIPT: Approximate MS(2)-AR(1) 
%}
clear all
close all
%% Model Parameters: y_t = \phi y_t-1 + \sig \varepsilon_t
phi_1 = 0.9; %AR coefficient
phi_2 = 0.1;
sig = 0.5; %std. dev of shock
N = 20; %Number of nodes
nS = 2; %Number of Regimes
Pi = [0.95,0.05;0.05,0.95];
barX = 0;



%% Construct f(\bar(y)_k | x)

%get nodes and weights via gauss-quadrature
[u,w] = hernodes(N);
xNodes = sig*sqrt(2)*u + barX;
phiCell = cell(2,1);
phiCell{1} = phi_1;
phiCell{2} = phi_2;
v{1} = barX;
v{2} = barX;
piMat = NaN(nS*N); %Markov matrix for discretized process
for ii=1:nS*N
    tMat = NaN(nS*N,1);
    ttMat = [(1:N),(1:N)];
    for jj=1:nS*N
        %tMat(jj) = normpdf(u(jj)-(phi*u(ii)))*(w(jj)./omega(u(jj)));
        tMat(jj) = normpdf(u(ttMat(jj))-(phiCell{floor(((ii-1)/N))+1}*u(ttMat(ii))))*...
            Pi(floor(((ii-1)/N))+1,floor(((jj-1)/N))+1);
        tMat(jj) = exp(-((u(ttMat(jj))-phiCell{floor(((ii-1)/N))+1}*u(ttMat(ii))-barX)^2 -(u(ttMat(jj))-barX)^2))...
            *w(ttMat(jj))/sqrt(pi)*Pi(floor(((ii-1)/N))+1,floor(((jj-1)/N))+1);
    end
    sX = sum(tMat);
    piMat(ii,:) = tMat'/sX;
end

nSim =300;
nTs = 300;
tsCell = cell(nSim,1);
sCell = cell(nSim,1);
aMat = NaN(nSim,4);
llMat = NaN(nSim,1);
for ii=1:nSim
    [tsCell{ii},sCell{ii}] = simMSMC(xNodes,piMat,0,nTs,N,nS);
end
Sig = [1,1];
v = [0,0];
P = Pi;
A = [0.5,0.5];
lb = [-Inf,-Inf,0,0];
ub = [Inf,Inf,1,1];
options = optimoptions(@fmincon,'Display','iter','UseParallel',1)
parpool(4)
for ii=1:nSim
    [aMat(ii,:),llMat(ii)] = fmincon(@(A) fn(A,tsCell{ii},Sig,v),...
        [0.9,0.1,0.95,0.95],[],[],[],[],lb,ub,[],options);
end

for ii=1:nSim
   if aMat(ii,1)<aMat(ii,2)
       tA = [aMat(ii,2),aMat(ii,1),aMat(ii,4),aMat(ii,3)];
       aMat(ii,:) = tA;
   end
end
%{
%% Check via simulation
nSim = 10000;
nTs = 1000;

tsCell = cell(nSim);
%Simulate nSim series of nTs
for ii=1:nSim
    tsCell{ii} = simMC(u,piMat,0,nTs,N);
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
mdl = arima('ARLags',1,'Constant',0,'Variance',sig,'AR',phi);
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
%}
function [ts,sMat] = simMSMC(u,piMat,initInd,nTs,N,nS)
    %u :: Nx1 set of nodes
    %Pi :: NXN Markov transition matrix
    %initInd :: {0,1,...,N} initial node
    
    S = NaN(nTs,1);
    sMat = NaN(nTs,1);
    ts = NaN(nTs,1);
    ttMat = [(1:N),(1:N)];
    if initInd==0
        q = getStatMarkov(piMat);
        tmp = rand;
        tNdx = 1;
        trig=0;
        while (tNdx<=nS*N)&&(trig==0)
           if tmp<=sum(q(1:tNdx))
               S(1) = tNdx;
               ts(1) = u(ttMat(tNdx));
               sMat(1) = floor((tNdx-1)/N)+1;
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
        while (tNdx<=nS*N)&&(trig==0)
           if tmp<=sum(piMat(S(tt-1),1:tNdx))
               S(tt) = tNdx;
               ts(tt) = u(ttMat(tNdx));
               sMat(tt) = floor((tNdx-1)/N)+1;
               trig = 1;
           else
               tNdx = tNdx + 1;
           end
        end
    end            
end

function [fnX] = fn(A,simTs,Sig,v)
    if A(1)>=A(2)
        fnX = -1*blhkLlMSIAH(simTs,Sig,A(1:2),v,[A(3),1-A(3);1-A(4),A(4)]...
            ,1);
    else
        fnX = -1*blhkLlMSIAH(simTs,Sig,[A(2),A(1)],v,...
            [A(4),1-A(4);1-A(3),A(3)],1);
    end
end

        
    
    
    


