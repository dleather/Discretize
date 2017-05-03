function [TF] = chkMssMsvar(PhiCell,Pi)
%{ 
    This functions follows Cho (2016, RED) to check is the MSVAR is
    stable in the mean square stability sense. This is equivalent to
    checking that the spectral radius of the matrix 
    (eye(n^2)-R) is less than unity, 

    s.t.

        R = [\pi_{1,1}kron(Phi(1),Phi(1)),...,\pi_{N,1}kron(Phi(1),Phi(1));
             \pi_{1,2}kron(Phi(2),Phi(2)),...,\pi_{N,2}kron(Phi(2),Phi(2));
                       .                                   .
                       .                                   .
                       .                                   .
             \pi_{N,1}kron(Phi(1),Phi(1)),...,\pi_{N,N}kron(Phi(N),Phi(N))]

        \pi_{i,j} = Pr( S_{t+1}=j | S_t = i ) 

    IN 

        PhiCell :: (Ns x 1) cell of (N x N) autoregrssive matrices
        Pi :: (Ns x Ns) matrix s.t Pi_{i,j} = Pr( S_{t+1} = j | S_t = i)
        
    OUT 

        TF = 1 if max(abs(eig(R)))<1
        TF = 0 o.w.
    %}

    Ns = size(PhiCell,1);
    N = size(PhiCell{1},1);
    [nn,mm] = size(Pi);
    if (nn~=Ns)||(mm~=Ns)
        error('Pi/PhiCell dimension mismatch')    
    end
    
    
    R = NaN(Ns*N^2);
    for ii=1:Ns
        for jj=1:Ns
            R(1+(ii-1)*(N^2):(N^2)*ii,1+(jj-1)*(N^2):(N^2)*jj) = ...
                Pi(jj,ii) .* kron(PhiCell{ii},PhiCell{ii});
        end
    end
    
    sr = max(abs(eig(R)));
    if sr < 1
        TF = 1;
    else
        TF = 0;
    end
    
    
    

    
    


           
end