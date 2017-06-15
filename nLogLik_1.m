
function [nLl] = nLogLik_1(data,params,smooth,K,M)
    %K=2, M=2, MSIAH
    % params is (K*M+2*((K^2)*M)+M^2 x 1) vector
    %   [v(:)-1:4; A(:)-5:12;. Sig(:) - 13:18 , P(:) - 19:20]
    v = reshape(params(1:K*M),K,M);
    A = reshape(params(K*M+1:K*M+K*K*M),K,K*M);
    P = r2pmat(params(end-((M*(M-1))-1):end));
    sig = cell(1,M);
    for ii=1:M
       %sig{ii} = sqrtm([(params(K*M+K^2*M+1+((K*(K+1))/2)*(ii-1)))^2,...
    %              params(K*M+K^2*M+((K*(K+1)/2))*ii);...
    %                  params(K*M+K^2*M+((K*(K+1)/2))*ii),...
    %              (params(K*M+K^2*M+2+((K*(K+1))/2)*(ii-1)))^2]);
       [~,sig{ii}] = sphere2cov(...
            params(K*M+K^2*M+(ii-1)*(K*(K+1)/2)+1:K*M+K^2*M+ii*(K*(K+1)/2)));
    end
    Sig = cell2mat(sig);
    if A(1,1)>A(1,3)
        A = [A(1:2,3:4),A(1:2,1:2)];
        v = [v(:,2),v(:,1)];
        Sig = [Sig(:,3:4),Sig(:,1:2)];
        P = [params(end),1-params(end);1-params(end-1),params(end-1)];
    end

    ll = blhkLlMSIAH(data,Sig,A,v,P,smooth);
    nLl=-ll;
end
