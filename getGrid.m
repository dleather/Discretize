function [Pr_mat_key,Pr_mat_key_pos,Pr_mat_intervals,zbar] = ...
    getGrid(A1bar,A2bar,SIGMAbar,N,m,method)
    
%{
    
    getGrid constructs the multivariate grid (pr_mat_kay), the position on the
    univariate grid (pr_mat_key_pos), the unadjusted intervals used in the Genz
    function (pr_mat_intervals), and the univariate grids (zbar). The code was
    taken from the code provided for Terry & Knotek (2011, Economic Letters)
    found on Stephen Terry's website.

    Model: x_{t+1} = A2bar x_{t} + A1bar +  \epsilon_{t+1}

    IN
        A1Bar :: (K x 1) :: intercept s.t. K is dimension of VAR
        A2Bar :: (K x K) :: reduced form autrogressive matrix
        SIGMAbar :: (K x K) :: covariance matrix of \epsilon_{t+1}
        N :: (K x 1) :: vector of number of discrete states per dimension
        m :: scalar :: number of standard deviations to form uniform grid
        method :: {0,1} :: method of grid formation
           1 -- Tauchen (1986) grid formation (uniform)
           2 -- Adda & Cooper grid formation (uni prob)

    OUT

        pr_mat_key :: (K x prod(N)) :: each column is a node of mvgrid
        pr_mat_key_pos :: (K x prod(N)) :: index of mv grid
        pr_mat_key_intervals :: (K x prod(N) x 2) :: unadjusted intervals for
            input into Genz (1992) code
        zbar :: (K x max(N)) :: matrix of univariate grids
    

  %}

    n = size(N,1);  %number of variables in VAR
    sstate_mean = (eye(n)-A2bar)\A1bar;
    SIGMAprocess = reshape((eye(n^2)-kron(A2bar,A2bar))\SIGMAbar(:),n,n);
    zbar = zeros(n,max(N));
    grid_stdev = diag(SIGMAprocess).^0.5;
    if method==1
        grid_increment = zeros(n,1);
        for i = 1:n
             grid_increment(i) = 2*m*grid_stdev(i)/(N(i)-1);
             zbar(i,1) = -m*grid_stdev(i) + sstate_mean(i);
            for j = 1:N(i)-1
                zbar(i,j+1) = zbar(i,j) + grid_increment(i);
            end
        end
    elseif method==2
        d = zeros(n,max(N));
        b = -4:.005:4;
        c = normcdf(b,0,1);
        for i = 1:n
            a = (1/(2*N(i))):(1/N(i)):1;
            for j = 1:N(i)
                [d1,d(i,j)] = min((a(j)-c).^2);
            end
            zbar(i,1:N(i)) = grid_stdev(i)*b(d(i,:))+sstate_mean(i);
        end
    end

    %compute key matrix & pos matrix
    Pr_mat_key = zeros(length(N),prod(N));
    Pr_mat_key_pos = zeros(length(N),prod(N));
    Pr_mat_key(length(N),:) = repmat(zbar(length(N),1:N(length(N))),...
        [1 prod(N)/N(length(N))]);
    Pr_mat_key_pos(length(N),:) = repmat(1:N(length(N)),[1 prod(N)/N(length(N))]);
    for i=length(N)-1:-1:1
        Pr_mat_key(i,:) = ...
            repmat(kron(zbar(i,1:N(i)),ones(1,prod(N(i+1:length(N))))),...
            [1 prod(N)/prod(N(i:length(N)))]);
        Pr_mat_key_pos(i,:) = ...
            repmat(kron(1:N(i),ones(1,prod(N(i+1:length(N))))),...
            [1 prod(N)/prod(N(i:length(N)))]);
    end

    nstate = prod(N);
    Pr_mat_intervals = zeros(n,nstate,2);   %this will store the unadjusted
        %limits of integration for each variable in each state, for input
        %into the Genz code
    if method==1
        for i = 1:nstate   %number of states
            for j = 1:n    %number of variables
                if Pr_mat_key_pos(j,i)==1
                    Pr_mat_intervals(j,i,1) = -inf;
                    Pr_mat_intervals(j,i,2) = zbar(j,Pr_mat_key_pos(j,i)) ...
                        + (grid_increment(j)/2);
                elseif Pr_mat_key_pos(j,i)==N(j)
                    Pr_mat_intervals(j,i,1) = zbar(j,Pr_mat_key_pos(j,i)) - ...
                        (grid_increment(j)/2);
                    Pr_mat_intervals(j,i,2) = inf;
                else
                    Pr_mat_intervals(j,i,1) = zbar(j,Pr_mat_key_pos(j,i)) - ...
                        (grid_increment(j)/2);
                    Pr_mat_intervals(j,i,2) = zbar(j,Pr_mat_key_pos(j,i)) + ...
                        (grid_increment(j)/2);
                end
            end
        end
    elseif method==2
        for i = 1:nstate  %number of states
            for j = 1:n    %number of variables
                if Pr_mat_key_pos(j,i)==1
                    Pr_mat_intervals(j,i,1) = -inf;
                    Pr_mat_intervals(j,i,2) = zbar(j,Pr_mat_key_pos(j,i)) + ...
                    (zbar(j,Pr_mat_key_pos(j,i)+1)-zbar(j,Pr_mat_key_pos(j,i)))/2;
                elseif Pr_mat_key_pos(j,i)==N(j)
                    Pr_mat_intervals(j,i,1) = zbar(j,Pr_mat_key_pos(j,i)) - ...
                    (zbar(j,Pr_mat_key_pos(j,i))-zbar(j,Pr_mat_key_pos(j,i)-1))/2;
                    Pr_mat_intervals(j,i,2) = inf;
                else
                    Pr_mat_intervals(j,i,1) = zbar(j,Pr_mat_key_pos(j,i)) - ...
                    (zbar(j,Pr_mat_key_pos(j,i))-zbar(j,Pr_mat_key_pos(j,i)-1))/2;
                    Pr_mat_intervals(j,i,2) = zbar(j,Pr_mat_key_pos(j,i)) + ...
                    (zbar(j,Pr_mat_key_pos(j,i)+1)-zbar(j,Pr_mat_key_pos(j,i)))/2;
                end
            end
        end
    end
