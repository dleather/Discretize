% Code to test approximations of AR(1) processes
% ... uses Tauchen's (Economic Letters, 1986) method
% ... and Tauchen and Hussey's (Econometrica, 1991) method
%
% Note: The tauchenHussey function uses Mario Miranda and Paul Fackler's
% function qnwnorm in their CompEcon toolbox.
% CompEcon is available at http://www4.ncsu.edu/~pfackler/compecon/
%
% For further details and results, see Flodén (2007), "A Note on the Accuracy 
% of Markov-Chain Approximations of Highly Persistent AR(1)-Processes", Stockholm
% School of Economics
%
% Martin Flodén, 2007
% martin.floden@hhs.se


clear all

% =========================================================================
% USER SETUP. Set PROCESS = 'XXX' to loop over many rho
PROCESS   = 'HSZ';              % Income process: AIY, HSZ, or STY, XXX
N         = 15;                 % # of nodes
GRAPHS    = 1;                  % graphs on (1) or off (0)
%METHODS                        % 1 = Tauchen
                                % 2 = Tauchen & Hussey, conditional variance
                                % 3 = Tauchen & Hussey, unconditional variance
                                % 4 = Tauchen & Hussey, weighted conditional and unconditional variance
                                % 5 = Adda & Cooper
% =========================================================================

% Income process
if PROCESS == 'AIY'
    loop_rho = 0.6000;
    sigma    = sqrt(0.013);
elseif PROCESS == 'HSZ'
    loop_rho = 0.9500;
    sigma    = sqrt(0.030);
elseif PROCESS == 'STY'
    loop_rho = 0.9800;
    sigma    = sqrt(0.020);
else
    loop_rho = [[0.5:0.1:0.9] [0.91:0.01:0.98]];
    sigma = 0.10;
end


for rho = loop_rho

    sigmaZ   = sigma/sqrt(1-rho^2);

    RHO    = zeros(5,1);
    CVZ    = zeros(5,1);
    UVZ    = zeros(5,1);
    CHAINS = [];


    for METHOD = 1:5

        if METHOD == 1
            [Z,PI] = tauchen(N,0,rho,sigma,1.2*log(N));
        elseif METHOD == 2
            [Z,PI] = tauchenhussey(N,0,rho,sigma,sigma);
        elseif METHOD == 3
            [Z,PI] = tauchenhussey(N,0,rho,sigma,sigmaZ);
        elseif METHOD == 4
            w      = 0.5 + rho/4;
            [Z,PI] = tauchenhussey(N,0,rho,sigma,w*sigma + (1-w)*sigmaZ);
        elseif METHOD == 5
            [Z,PI] = addacooper(N,0,rho,sigma);
        else
            error('Method not implemented')
        end

        % Calculate ergodic distribution with method suggested by Paul Klein
        % ... just using PI^100000 etc will not work for highly persistent processes
        [evec,eval] = eig(PI');
        [tmp,i]     = min(abs(1-diag(eval)));
        evec        = evec(:,i);
        ergPI       = evec / sum(evec);

        CHAINS(METHOD).Z     = Z;
        CHAINS(METHOD).PI    = PI;
        CHAINS(METHOD).ergPI = ergPI;


        Zprime = ones(N,1)*Z';
        dz     = Z*ones(1,N) - Zprime;

        EZ          = ergPI'*Z;                                     % expected Z
        UVZ(METHOD) = ergPI'*((Z-EZ).^2);                           % unconditional variance of Z
        CEZ         = PI*Z;                                         % conditional expectation of Z
        vecCVZ      = sum(PI.*((Zprime-CEZ*ones(1,N)).^2),2);       % conditional variance in each node
        CVZ(METHOD) = ergPI'*vecCVZ;                                % weighted conditional variance

        zpz         = (Zprime-EZ).*(Z*ones(1,N)-EZ);
        vecRHO      = sum(PI.*zpz,2)./((Z-EZ).^2);
        ii          = find(abs(Z-EZ)>1e-10);
        jj          = find(abs(Z-EZ)<1e-10);
        RHO(METHOD) = ergPI(ii)'*vecRHO(ii)/sum(ergPI(ii));
        vecRHO(jj)  = NaN;

        if GRAPHS
            figure
            subplot(1,2,1)
            plot(Z,vecRHO,'-x',Z,rho*ones(size(Z)),'r')
            xlabel('Z')
            ylabel('\rho')
            title('\rho conditional on node (blue) and true (red)')

            subplot(1,2,2)
            plot(Z,sqrt(vecCVZ),'-x',Z,sigma*ones(size(Z)),'r')
            xlabel('Z')
            ylabel('\sigma_\epsilon')
            title('\sigma_\epsilon conditional on node (blue) and true (red)')
        end

    end

    fprintf('============================================================================\n')
    fprintf('Parameter                         Implied\n')
    fprintf('Parameter      True      Tauchen    TH(unc)   TH(cond)     TH(w)      AC\n')
    fprintf('rho         %9.4f  %9.4f  %9.4f  %9.4f  %9.4f  %9.4f\n',rho,RHO)
    fprintf('std(eps)    %9.4f  %9.4f  %9.4f  %9.4f  %9.4f  %9.4f\n',sigma,sqrt(CVZ))
    fprintf('std(z)      %9.4f  %9.4f  %9.4f  %9.4f  %9.4f  %9.4f\n',sigmaZ,sqrt(UVZ))
    fprintf('----------------------------------------------------------------------------\n')
    fprintf('Parameter                      generated/true\n')
    fprintf('Parameter                Tauchen    TH(unc)   TH(cond)     TH(w)      AC\n')
    fprintf('rho                    %9.4f  %9.4f  %9.4f  %9.4f  %9.4f\n',RHO/rho)
    fprintf('std(eps)               %9.4f  %9.4f  %9.4f  %9.4f  %9.4f\n',sqrt(CVZ)/sigma)
    fprintf('std(z)                 %9.4f  %9.4f  %9.4f  %9.4f  %9.4f\n',sqrt(UVZ)/sigmaZ)
    fprintf('----------------------------------------------------------------------------\n')
    fprintf('max(node)/std(z)     ')
    for i = 1:5
        fprintf('  %9.4f',CHAINS(i).Z(end)/sigmaZ)
    end
    fprintf('\nmin(erg)/max(erg)    ')
    for i = 1:5
        fprintf('  %9.4f',min(CHAINS(i).ergPI)/max(CHAINS(i).ergPI))
    end
    fprintf('\n');

end