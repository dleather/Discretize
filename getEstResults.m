%{
    getEstResults.m

    This script should create a table and figures for the estimation
%}  

load('assetPriceResults_cluster.mat')

%% Create Table

%{ Params | true value | mean w/ 95% CI underneathe | median }

% Get matrix of means

estMat = cell2mat(estCell');
dEstMat = cell2mat(dEstCell');
meanEst = mean(estMat,2);
dMeanEst = mean(dEstMat,2);
medEst = median(estMat,2);
dMedEst = median(dEstMat,2);
lbPrcEst = prctile(estMat,2.5,2);
ubPrcEst = prctile(estMat,97.5,2);
dLbPrcEst = prctile(dEstMat,2.5,2);
dUbPrcEst = prctile(dEstMat,97.5,2);

dataMat = [meanEst, medEst, lbPrcEst, ubPrcEst, dMeanEst, dMedEst, dLbPrcEst, dUbPrcEst];

input.tableColumnAlignment = 'c';
input.tableCaption = 'Mean, Median, and 95\% Percentiles of 120 simulations from a continuous and discretized MS(2)-VAR(1), $x_t = \Phi(S_t) x_{t-1} + \mu(S_t)+\varepsilon_t$ s.t. $\varepsilon_t \sim N(0,\Sigma(S_t))$';
input.tableLabel = 'Discretization Test';
input.tableRowLabels = {'$\mu(1)_1$','$\mu(1)_2$','$\mu(2)_1$',...
    '$\mu(2)_2$','$\Phi(1)_{1,1}$','$\Phi(1)_{2,1}$','$\Phi(1)_{1,2}$'...
    '$\Phi(1)_{2,2}$','$\Phi(2)_{1,1}$','$\Phi(2)_{2,1}$','$\Phi(2)_{1,2}$'...
    '$\Phi(2)_{2,2}$','$\sigma(1)^2_1$','$\sigma(1)^2_2$','$\rho(1)_{1,2}$',...
    '$\sigma(2)^2_1$','$\sigma(2)^2_2$','$\rho(2)_{1,2}$','$\pi_{1,1}$','$\pi_{2,2}$'};
input.tableColLabels = {'Mean','Median','$2.5^{th}$ Percentile',...
    '$97.5^{th}$ Percentile','Mean','Median','$2.5^{th}$ Percentile',...
    '$97.5^{th}$ Percentile'};
input.data = dataMat;
latex = latexTable(input);





%% Create Figures

% For each main group, has a panel of overlaying histograms w/ true value
% as vline
% Groups: AR Coefficients (8), Constant (4), VarCov (6), Probability (2)

%Group 1: Constants (4)
figure('Units','inches',...
'Position',[0 0 7.5 10],...
'PaperPositionMode','auto');
subplot(2,2,1);
histogram(estMat(1,:),'BinMethod','fd','Normalization','probability')
hold on
histogram(dEstMat(1,:),'BinMethod','fd','Normalization','probability');
vline(muCell{1}(1),'-+k')
legend('Continuous','Discretized')
title('$\mu(1)_1$','Interpreter','latex')
subplot(2,2,3);
histogram(estMat(2,:),'BinMethod','fd','Normalization','probability')
hold on
histogram(dEstMat(2,:),'BinMethod','fd','Normalization','probability');
vline(muCell{1}(2),'-+k')
title('$\mu(1)_2$','Interpreter','latex')
subplot(2,2,2);
histogram(estMat(3,:),'BinMethod','fd','Normalization','probability')
hold on
histogram(dEstMat(3,:),'BinMethod','fd','Normalization','probability');
vline(muCell{2}(1),'-+k')
title('$\mu(2)_1$','Interpreter','latex')
subplot(2,2,4);
histogram(estMat(4,:),'BinMethod','fd','Normalization','probability')
hold on
histogram(dEstMat(4,:),'BinMethod','fd','Normalization','probability');
vline(muCell{2}(2),'-+k')
title('$\mu(2)_2$','Interpreter','latex')
suptitle('Approximation of Estimate Distrbution, \mu(S_t)')
print -depsc2 fig1.eps

%Group 2: Ar Coefficients
figure('Units','inches',...
'Position',[0 0 10 7.5],...
'PaperPositionMode','auto');
subplot(2,4,1)
histogram(estMat(5,:),'BinMethod','fd','Normalization','probability')
hold on
histogram(dEstMat(5,:),'BinMethod','fd','Normalization','probability')
vline(PhiCell{1}(1,1),'-+k')
title('$\Phi(1)_{1,1}$','Interpreter','latex')
legend('Continuous','Discrete')
subplot(2,4,5)
histogram(estMat(6,:),'BinMethod','fd','Normalization','probability')
hold on
histogram(dEstMat(6,:),'BinMethod','fd','Normalization','probability')
vline(PhiCell{1}(2,1),'-+k')
title('$\Phi(1)_{2,1}$','Interpreter','latex')
subplot(2,4,2)
histogram(estMat(7,:),'BinMethod','fd','Normalization','probability')
hold on
histogram(dEstMat(7,:),'BinMethod','fd','Normalization','probability')
vline(PhiCell{1}(1,2),'-+k')
title('$\Phi(1)_{1,2}$','Interpreter','latex')
subplot(2,4,6)
histogram(estMat(8,:),'BinMethod','fd','Normalization','probability')
hold on
histogram(dEstMat(8,:),'BinMethod','fd','Normalization','probability')
vline(PhiCell{1}(2,2),'-+k')
title('$\Phi(1)_{2,2}$','Interpreter','latex')
subplot(2,4,3)
histogram(estMat(9,:),'BinMethod','fd','Normalization','probability')
hold on
histogram(dEstMat(9,:),'BinMethod','fd','Normalization','probability')
vline(PhiCell{2}(1,1),'-+k')
title('$\Phi(2)_{1,1}$','Interpreter','latex')
subplot(2,4,7)
histogram(estMat(10,:),'BinMethod','fd','Normalization','probability')
hold on
histogram(dEstMat(10,:),'BinMethod','fd','Normalization','probability')
vline(PhiCell{2}(2,1),'-+k')
title('$\Phi(2)_{2,1}$','Interpreter','latex')
subplot(2,4,4)
histogram(estMat(11,:),'BinMethod','fd','Normalization','probability')
hold on
histogram(dEstMat(11,:),'BinMethod','fd','Normalization','probability')
vline(PhiCell{2}(1,2),'-+k')
title('$\Phi(2)_{1,2}$','Interpreter','latex')
subplot(2,4,8)
histogram(estMat(12,:),'BinMethod','fd','Normalization','probability')
hold on
histogram(dEstMat(12,:),'BinMethod','fd','Normalization','probability')
vline(PhiCell{2}(2,2),'-+k')
title('$\Phi(2)_{2,2}$','Interpreter','latex')
suptitle('Approximation of Estimate Distrbution, \Phi(S_t)')
print -depsc2 fig2.eps

% Group 3: Covariance Matrices

%Translate from sphere2cov
for ii=1:size(estMat,2)
   tmp = sphere2cov(estMat(13:15,ii));
   estMat(13,ii) = tmp(1,1);
   estMat(14,ii) = tmp(2,2);
   estMat(15,ii) = tmp(1,2);
   tmp = sphere2cov(dEstMat(13:15,ii));
   dEstMat(13,ii) = tmp(1,1);
   dEstMat(14,ii) = tmp(2,2);
   dEstMat(15,ii) = tmp(1,2);
   tmp = sphere2cov(estMat(16:18,ii));
   estMat(16,ii) = tmp(1,1);
   estMat(17,ii) = tmp(2,2);
   estMat(18,ii) = tmp(1,2);
   tmp = sphere2cov(dEstMat(16:18,ii));
   dEstMat(16,ii) = tmp(1,1);
   dEstMat(17,ii) = tmp(2,2);
   dEstMat(18,ii) = tmp(1,2);
end

%get covariance matrices
sigCell{1} = sigCell{1}*sigCell{1}';
sigCell{2} = sigCell{2}*sigCell{2}';

figure('Units','inches',...
'Position',[0 0 7.5 10],...
'PaperPositionMode','auto');
subplot(3,2,1)
histogram(estMat(13,:),'BinMethod','fd','Normalization','probability')
hold on
histogram(dEstMat(13,:),'BinMethod','fd','Normalization','probability')
vline(sigCell{1}(1,1),'-+k')
title('$\Sigma(1)_{1,1}$','Interpreter','latex')
legend('Continuous','Discrete')
subplot(3,2,3)
histogram(estMat(14,:),'BinMethod','fd','Normalization','probability')
hold on
histogram(dEstMat(14,:),'BinMethod','fd','Normalization','probability')
vline(sigCell{1}(2,2),'-+k')
title('$\Sigma(1)_{2,2}$','Interpreter','latex')
subplot(3,2,5)
histogram(estMat(15,:),'BinMethod','fd','Normalization','probability')
hold on
histogram(dEstMat(15,:),'BinMethod','fd','Normalization','probability')
vline(sigCell{1}(1,2),'-+k')
title('$\Sigma(1)_{1,2}$','Interpreter','latex')
subplot(3,2,2)
histogram(estMat(16,:),'BinMethod','fd','Normalization','probability')
hold on
histogram(dEstMat(16,:),'BinMethod','fd','Normalization','probability')
vline(sigCell{2}(1,1),'-+k')
title('$\Sigma(2)_{1,1}$','Interpreter','latex')
subplot(3,2,4)
histogram(estMat(17,:),'BinMethod','fd','Normalization','probability')
hold on
histogram(dEstMat(17,:),'BinMethod','fd','Normalization','probability')
vline(sigCell{2}(2,2),'-+k')
title('$\Sigma(2)_{2,2}$','Interpreter','latex')
subplot(3,2,6)
histogram(estMat(18,:),'BinMethod','fd','Normalization','probability')
hold on
histogram(dEstMat(18,:),'BinMethod','fd','Normalization','probability')
vline(sigCell{2}(1,2),'-+k')
title('$\Sigma(2)_{1,2}$','Interpreter','latex')
suptitle('Approximation of Estimate Distrbution, \Sigma(S_t)')
print -depsc2 fig3.eps

%Group 4: Probabilties
estMat(19,:) = 1 ./ (1 + exp(estMat(19,:)));
estMat(20,:) = 1 ./ (1 + exp(estMat(20,:)));
dEstMat(19,:) = 1 ./ (1 + exp(dEstMat(19,:)));
dEstMat(20,:) = 1 ./ (1 + exp(dEstMat(20,:)));
figure('Units','inches',...
'Position',[0 0 5 7.5],...
'PaperPositionMode','auto');
subplot(2,1,1)
histogram(estMat(19,:),'BinMethod','fd','Normalization','probability')
hold on
histogram(dEstMat(19,:),'BinMethod','fd','Normalization','probability')
vline(Pi(1,1),'-+k')
title('$\pi_{1,1}$','Interpreter','latex')
legend('Continuous','Discrete')
subplot(2,1,2)
histogram(estMat(20,:),'BinMethod','fd','Normalization','probability')
hold on
histogram(dEstMat(20,:),'BinMethod','fd','Normalization','probability')
vline(Pi(2,2),'-+k')
title('$\pi_{2,2}$','Interpreter','latex')
suptitle('Approximation of Estimate Distribution,\Pi')
print -depsc2 fig4.eps






