function Pe=BSC(N,Eu)
addpath('D:\studing\yanjiusheng\后处理\post-processing-m\2IR\Encode\BSC')
%%% LL June 11, 2015
% calculate the error probability for each bit channel
% sort the bit channels according to their error probability
% select the best bit channels as the information set
n=log2(N);
v = 2^6;
u = 2*v; 

W = [1-Eu, Eu]; % W stores W(y|0) W(\bar(y)|0)
   
%%
Pe=zeros(1,N);
% given a bit channel i
parfor i=1:N
    fprintf('\nNow iter : %3d',i);
    % Initial processing 
    W1 = sortTran_sim(W);
    [Q] = degrading_merge(W1, u);
    bi = dec2bin(i-1,n);
    for m=1:n 
        % sort the LRs according to bi(m)
        if bi(m) == '0'
            % size of the output alphabet at this level  
            [W10] = calcTran_sim_combine(Q,0); % type 0 transformation 
            % calculate LRs for W
            [W20] = sortTran_sim(W10);
            [Q] = degrading_merge(W20,u);
        else
            [W11] = calcTran_sim_combine(Q,1); % type 1 transformation 
            % calculate LRs for W
            [W21] = sortTran_sim(W11);
            [Q] = degrading_merge(W21,u);
        end
    end
    nq = length(Q);
    Pe(i) = sum(Q(nq/2+1:nq));
end
name1='%d_%G.mat';
name=sprintf(name1,N,Eu);
save(name,'Pe');
end