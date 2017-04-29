function [ts,S] = simMC(u,Pi,initInd,nTs,N)
    %u :: NxK set of nodes
    %Pi :: NXN Markov transition matrix
    %initInd :: {0,1,...,N} initial node
    
    S = NaN(nTs,1);
    ts = NaN(nTs,size(u,2));
    if initInd==0
        q = getStatMarkov(Pi);
        tmp = rand;
        tNdx = 1;
        trig=0;
        while (tNdx<=N)&&(trig==0)
           if tmp<=sum(q(1:tNdx))
               S(1) = tNdx;
               ts(1,:) = u(tNdx,:);
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
               ts(tt,:) = u(tNdx,:);
               trig = 1;
           else
               tNdx = tNdx + 1;
           end
        end
    end            
end


        