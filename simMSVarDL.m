%{
simMSVarDL
    Author: David Leather -- Economics Department -- UNC Chapel Hill
    funtion: Simulate a MS-VAR(1) w/ iid Gaussian errors
    IN
        muCell :: (nRegime x 1) cell w/ constant vector
        PhiCell :: (nRegime x 1) cell w/ autoregressive matrix
        sigCell :: (nRegime x 1) cell w/ sqrt of Covariance matrix
        Pi :: Transition matrix
        nTs :: Length of time series to be simulated
        nSim :: Numbers of simulations
        q :: vector that represents distribution of original state
        init :: Vector of start values
    OUT
        simCell :: Cell of simulations
        sCell :: Cell w/ simualted hidden Markov regime paths
%}

function [simCell,sCell] = simMSVarDL(muCell,PhiCell,sigCell,nTs,nSim,Pi,q,init,initX)
%% Input Check
    %set.seed = seed;
    nRegime = size(muCell,1);
    if (size(PhiCell,1)~=nRegime)
        error('Size of PhiCell & muCell not the same')
    elseif (size(sigCell,1)~=nRegime)
        error('Size of SigCell & MuCell different')
    elseif ~(isscalar(nTs))
        error('nTs is not scalar')
    elseif ~(isscalar(nSim))
        error('nSim is not scalar')
    elseif max(size(q,1),size(q,2))~=nRegime
        error('muCell dimensions and q vector dimensions don''t match up')
    end
    
    nVar = size(PhiCell{1},1);
    
%% Define Structures
    sCell = cell(nSim,1);
    simCell = cell(nSim,1);
    epsCell = cell(nSim,1);
    
    for iSim = 1:nSim
        simCell{iSim} = NaN(nVar,nTs+1);
        sCell{iSim} = NaN(nTs+1,1);
        epsCell{iSim} = randn(nVar,nTs);
        
%% Draw initial regime & Set initial data
        tRand = rand(1,1);
        temp = 0;
        tSwitch = 0;
        iR = 1;
		if init==0
			while((tSwitch==0)||(iR~=nRegime+1))
			   temp = temp + q(iR);
			   if tRand<=temp
				   sCell{iSim}(1,1) = iR;
				   tSwitch = 1;
			   end
			   iR = iR + 1;
			end
		else
			sCell{iSim}(1,1) = init;
        end
        
        simCell{iSim}(:,1) = initX;
      
%%  Simulate series for each simulation
        for t=1:nTs
           %Draw new state
           tRand = rand(1,1);
            temp = 0;
            tSwitch = 0;
            iR = 1;
            while((tSwitch==0)&&(iR~=nRegime+1))
               temp = temp + Pi(sCell{iSim}(t),iR);
               if tRand<=temp
                   sCell{iSim}(t+1) = iR;
                   tSwitch = 1;
               end
               iR = iR + 1;
            end
            
            %Generate next sequence
            simCell{iSim}(:,t+1) = muCell{sCell{iSim}(t+1)} + ...
                PhiCell{sCell{iSim}(t+1)}*simCell{iSim}(:,t) + ...
                sigCell{sCell{iSim}(t+1)}*epsCell{iSim}(:,t);
        end
end