function [Pr_mat] = ...
    getPrMat(A1bar,A2bar,SIGMAbar,Pr_mat_intervals_j,Pr_mat_key_i,...
        nstate_i,nstate_j)
    n = size(A1bar,1);
    error_est = zeros(nstate_i,nstate_j);
    Pr_mat_intervals_adjusted = zeros(n,nstate_j,2);
    Pr_mat = zeros(nstate_i,nstate_j);

    random_draws = 300;

    for i = 1:nstate_i; %rows of Pr_mat
        Pr_mat_intervals_adjusted(:,:,1) = Pr_mat_intervals_j(:,:,1) - ...
            repmat((A1bar + A2bar*Pr_mat_key_i(:,i)),1,nstate_j);
        Pr_mat_intervals_adjusted(:,:,2) = Pr_mat_intervals_j(:,:,2) - ...
            repmat((A1bar + A2bar*Pr_mat_key_i(:,i)),1,nstate_j);
        for j = 1:nstate_j;   %columns of Pr_mat
            %Pr_mat(i,j) = P(state j|state i)
            [Pr_mat(i,j), error_est(i,j)] = ...
                qscmvnv(random_draws,SIGMAbar,Pr_mat_intervals_adjusted(:,j,1),...
                eye(n),Pr_mat_intervals_adjusted(:,j,2));
        end;
    end;

    %rounding error adjustment
    round_sum = sum(Pr_mat,2);
    for i = 1:size(Pr_mat,1);
        Pr_mat(i,:) = Pr_mat(i,:)/round_sum(i);
    end
end
