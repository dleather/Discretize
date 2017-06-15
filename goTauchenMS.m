%{ 
SCRIPT :: goTauchenMS.m 

This code should accomplish the following :
    
    1) Discretize each specification of MS-VAR
    2) Run Tests
    2) Simultate nSim discrete processes
    3) Simulate nSim continuous processes
    4) Compare means, median and confidence interval of estimates
   
TODO

    1) Test Asset Pricing capabilities

PARAMETERS

    nSim :: number of simulations of discrete & continuous DGP
    nTs :: length of each simulated series.
    N :: (K x Ns) matrix of nummber of discretized nodes for each
        dimension, and regime
    PhiCell :: (Ns x 1) cell of (K x K) autoregressive matrices
    muCell :: (Ns x 1) cell of (K x 1) intercept matrices
    CovCell :: (Ns x 1) cell of (K x K) covariance matrices
    Pi :: (Ns x Ns) transition matrix s.t. rows sum to unity
    m :: Number of std. dev. from mean to form grid
    
SPECIFICATIONS
    K = 2,3
    Ns = 2,4,8
    N = anisotropic grid, isotropic grid

%}
clear all
close all

simACell = cell(15,1);
cmpACell = cell(15,1);
errCell = cell(15,1);
ticMat = NaN(15,1);
tmpRandCell = cell(15,1);
for iN=15:15
%% Global parameters
nSim = 10;
nTs = 10000;
m = 3;
llCell = cell(nSim,1);
nnSim = 100;
nnTs = 1200;

for method = 1:1
   %% Set local parameters
   if method==1
        K = 2;
        Ns = 2;
        N = [iN+4,iN+4;iN+4,iN+4];
        PhiCell = cell(Ns,1);
        PhiCell{1} = [0.5, 0.1; 0, 0.7];
        PhiCell{2} = [0.85, 0; 0.1 , 0.2];
        muCell = cell(Ns,1);
        muCell{1} = [0.3; 0.2];
        muCell{2} = [0.8; 00.5];
        CovCell = cell(Ns,1);
        CovCell{1} = eye(K);
        CovCell{2} = [0.5, 0.2; 0.2, 0.8];
        Pi = [ 0.95, 0.05; 0.09, 0.91 ];
   else
       error('you''ve ran out of methods')
   end
   
   %Check for MSS
   MSS = chkMssMsvar(PhiCell,Pi);
   if MSS==0
       error('MSVAR not MSS')
   end
   
   %% Discretization
   tic;
   prMatKeyCell = cell(Ns,1);
   prMatKeyPosCell = cell(Ns,1);
   prMatIntCell = cell(Ns,1);
   zBarCell = cell(Ns,1);

   % Get multivariate grids
   for jj=1:Ns
       [prMatKeyCell{jj},prMatKeyPosCell{jj},prMatIntCell{jj},zBarCell{jj}]...
           = getGrid(muCell{jj},PhiCell{jj},CovCell{jj},N(:,jj),m,1);
   end
   
   % Get Pi_{i,j}
   PiIJCell = cell(Ns,1);
   for ii=1:Ns
       for jj=1:Ns
           PiIJCell{ii,jj} = ...
               getPrMat(muCell{jj},PhiCell{jj},CovCell{jj},prMatIntCell{jj},...
                    prMatKeyCell{ii},size(prMatKeyCell{ii},2),...
                    size(prMatKeyCell{jj},2));
       end
   end
   
   % Mix Pi_{i,j}s to form dPi
   NNVec = NaN(Ns,1);  %Total number of discrete states between regimes
   for ii=1:Ns
       NNVec(ii) = size(prMatKeyCell{ii},2);
   end
   cumNN = cumsum(NNVec);
   NN = sum(NNVec);
   ndx = NaN(NN,1); %This vector indexes that state at each discrete state
   for ii=1:Ns
       if ii~=1
            ndx(1+cumNN(ii-1):cumNN(ii)) = ii;
       else
           ndx(1:cumNN(ii)) = ii;
       end
   end
   
   dPi = NaN(NN);
   cnt = 1;
   for ii=1:NN
       for jj=1:Ns
           if jj~=1
                dPi(ii,1+cumNN(jj-1):cumNN(jj)) = ...
                    Pi(ndx(ii),jj).*PiIJCell{ndx(ii),jj}(cnt,:);
           else
               dPi(ii,1:NNVec(jj)) = ...
                   Pi(ndx(ii),jj).*PiIJCell{ndx(ii),jj}(cnt,:);
           end
       end
       if ii~=NN
           if ndx(ii)==ndx(ii+1)
               cnt = cnt + 1;
           else
               cnt = 1;
           end
       end
   end
   ticMat(iN) = toc;
   %% Simulate nSim of each discrete & continuous series series
   dTsCell = cell(nSim,1);
   dSCell = cell(nSim,1);
   tsCell = cell(nSim,1);
   sCell = cell(nSim,1);
   u = [];
   for ii=1:Ns
      u = [u;prMatKeyCell{ii}']; 
   end
   
   %Discrete simulation
   for ii=1:nSim
        [dTsCell{ii},dSCell{ii}] = simMC(u,dPi,0,nTs,NN);
        disp(ii)
   end
   
   %Continuous simulation
   sigCell = cell(Ns,1);
   for ii=1:Ns
      sigCell{ii} = sqrtm(CovCell{ii});
   end
   qq = getStatMarkov(Pi);
   [tsCell,sCell{ii}] = ...
        simMSVarDL(muCell,PhiCell,sigCell,10*nTs,nSim,Pi,qq,0,u(1,:)');
  
  %Estimate via MLE
   init = [muCell{1};muCell{2};PhiCell{1}(:);PhiCell{2}(:);sqrt(diag(sigCell{1}));...
           sigCell{1}(2,1);sqrt(diag(sigCell{2}));sigCell{2}(2,1);diag(Pi)];
   options = optimoptions(@fminunc,'Display','iter','StepTolerance',1e-06);
   dEstCell = cell(nSim,1);
   dLlMat = NaN(nSim,1);
   estCell = cell(nSim,1);
   llMat = NaN(nSim,1);
   for ii=1:nSim
        [dEstCell{ii}, dLlMat(ii)] = ...
            fminunc(@(X) nLogLik_1(dTsCell{ii},X,1,2,2),init,options);
        if dEstCell{ii}(5)>dEstCell{ii}(9)
            dEstCell{ii} = [dEstCell{ii}(3:4);dEstCell{ii}(1:2);...
                            dEstCell{ii}(9:12);dEstCell{ii}(5:8);...
                            dEstCell{ii}(16:18);dEstCell{ii}(13:15);...
                            dEstCell{ii}(20);dEstCell{ii}(19)];
        end
        [estCell{ii}, llMat(ii)] = ...
            fminunc(@(X) nLogLik_1(tsCell{ii}(:,end-(nTs+1):end)',X,1,2,2),init,options);
        if estCell{ii}(5)>estCell{ii}(9)
            estCell{ii} = [estCell{ii}(3:4);estCell{ii}(1:2);...
                            estCell{ii}(9:12);estCell{ii}(5:8);...
                            estCell{ii}(16:18);estCell{ii}(13:15);...
                            estCell{ii}(20);estCell{ii}(19)];
        end
         
   end
   
   
   resMat = NaN(20,4); %Mean, median, 2.5 prctle, 97.5 prctle
   dResMat = NaN(20,4);
   resMat = [mean(cell2mat(estCell'),2),prctile(cell2mat(estCell'),50,2),...
             prctile(cell2mat(estCell'),2.5,2),...
             prctile(cell2mat(estCell'),95,2)];
   dResMat = [mean(cell2mat(dEstCell'),2),prctile(cell2mat(dEstCell'),50,2),...
             prctile(cell2mat(dEstCell'),2.5,2),...
             prctile(cell2mat(dEstCell'),95,2)];
         
    mseVec = NaN(2,2);
    mseVec(1,1) = sum((mean(cell2mat(dEstCell'),2)-init).^2)./(nSim-1);
    mseVec(2,1) = sum((prctile(cell2mat(dEstCell'),50,2)-init).^2)./(nSim-1);
    mseVec(1,2) = sum((mean(cell2mat(estCell'),2)-init).^2)./(nSim-1);
    mseVec(2,2) = sum((prctile(cell2mat(estCell'),50,2)-init).^2)./(nSim-1);
   
   %}
   %Asset pricing
   a = 0.9; %arbitrary
   b = 0.2;
   %c = 0; for
   te = [a,b];
   d = 0.003;
   eta = NaN(NN,1);
   cnt=1;
   for mm=1:Ns
       for ii=1:NNVec(mm)
           eta(cnt) = exp(-te*prMatKeyCell{mm}(:,ii) - d);
           cnt = cnt + 1;
       end
   end
   
   dAp = (eye(NN)-eta.*dPi)\eta;
   
   %simulate asset prices:
   %Choose nA nodes at random from multivariate grid. Simulate nsim series
   % of length 300, construct product of eta for each t, and then average
   % to get prices. Compare to the corresponding discretized prices.
   nA = 30;
   tmpRand = randperm(size(dPi,1),nA);
   mvGrid = [prMatKeyCell{1},prMatKeyCell{2}];
   mvGridA = mvGrid(:,tmpRand);
   cmpA = dAp(tmpRand);
   simA = NaN(nA,1);
   for ii=1:nA
      %Simulate nsim=1000 series of  length 300 at point A
      [tsCellA,sCellA] = ...
        simMSVarDL(muCell,PhiCell,sigCell,nnTs,nnSim,Pi,qq,ndx(tmpRand(ii)),u(tmpRand(ii),:)');
      
      %Get Vector of Eta
      etaMat = NaN(nnSim,nnTs);
      for jj=1:nnSim
          etaMat(jj,:) = -te*tsCellA{jj}(:,1:end-1)-d;
      end
      etaMat = exp(etaMat);
      if any(any(isnan(etaMat)))
          keyboard;
      end
      prodEta = cumprod(etaMat,2);
      simA(ii) = mean(sum(prodEta,2));
      disp(ii)
   end
   simACell{iN} = simA;
   cmpACell{iN} = cmpA;
   tmpRandCell{iN} = tmpRand;
   errCell{iN} = simA - cmpA;
end
%
 
      
    
   end
   
   save('assetPriceResults.mat');
   
   
       
   
    
 
