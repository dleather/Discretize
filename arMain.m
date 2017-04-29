%{
SCRIPT: Approximate AR(1) 
%}
clear all

%% Model Parameters: y_t = \phi y_t-1 + \sig \varepsilon_t
phi = 0.95; %AR coefficient
sig = 0.03; %std. dev of shock
barX = 0;
N = 5; %Number of nodes
%omega = @(x) normpdf(x); % Set omega to std. normal PDF

%% Construct f(\bar(y)_k | x)

%get nodes and weights via gauss-quadrature
[u,w] = hernodes(N);
xNodes = sqrt(2)*sig*u + barX; %GH Nodes based on \sig_{\epsilon}
%xNodes = sqrt(2)*(sig/sqrt(1-phi^2))*u + barX; %Use conditional variance
%sigW = 0.5 + phi/4;
%xNodes = sqrt(2)*(sigW*sig+(1-sigW)*(sig/sqrt(1-phi^2)))*u + barX;
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
nSim = 350;
nTs = 100000;

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




        
    
    
    


