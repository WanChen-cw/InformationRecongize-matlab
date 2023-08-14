function polar_info_esti_all = Copy_of_SCL_decoder(Ksb,N, L, Sd,Eu)%2018.1.7.14:16 Yu Y. R.
%% 输入分段
M=length(Ksb)/N;
llr(Ksb==0)=log2((1-Eu)/Eu);
llr(Ksb==1)=log2(Eu/(1-Eu));
llr=reshape(llr,N,M)';
%% polar译码前置
llr = bitrevorder(llr')';
lambda_offset = 2.^(0 : log2(N));
llr_layer_vec = get_llr_layer(N);
bit_layer_vec = get_bit_layer(N);
frozen_bits=(Sd(1,:)~=-1);%冻结bit位置
K=sum(Sd(1,:)==-1);%信息比特数量
m = log2(N);
polar_info_esti_all=Sd;
%%
%memory declared
lazy_copy = zeros(M,m, L);
%If Lazy Copy is used, there is no data-copy in the decoding process. 
%We only need to record where the data come from. Here,data refer to LLRs and partial sums.
%initialize
lazy_copy(:,:, 1) = 1;
%Lazy Copy is a relatively sophisticated operation for new learners of polar codes. 
%If you do not understand such operation, you can directly copy data.
%If you can understand lazy copy and you just start learning polar codes
%for just fews days, you are very clever,

P = zeros(M,N - 1, L); %Channel llr is public-used, so N - 1 is enough.
C = zeros(M,N - 1, 2 * L);%I do not esitimate (x1, x2, ... , xN), so N - 1 is enough.
u = zeros(M,K, L);%unfrozen bits that polar codes carry, including crc bits.
PM = zeros(M,L);%Path metrics

activepath = zeros(M,L);%Indicate if the path is active. '1'→active; '0' otherwise.
activepath(:,1) = 1;
cnt_u = 1;%information bit counter 
%decoding starts
%default: in the case of path clone, the origianl path always corresponds to bit 0, while the new path bit 1.
for phi = 0 : N - 1
    layer = llr_layer_vec(phi + 1);
    phi_mod_2 = mod(phi, 2);
    
    %% P updating
    for M_index=1:M
    for l_index = 1 : L
        if activepath(M_index,l_index) == 0
            continue;
        end
        switch phi%Decoding bits u_0 and u_N/2 needs channel LLR, so the decoding of them is separated from other bits. 
            case 0
                index_1 = lambda_offset(m);
                for beta = 0 : index_1 - 1
                    P(M_index,beta + index_1, l_index) = sign(llr(M_index,beta + 1)) * sign(llr(M_index,beta + index_1 + 1)) * min(abs(llr(M_index,beta + 1)), abs(llr(M_index,beta + index_1 + 1)));
                end
                for i_layer = m - 2 : -1 : 0
                    index_1 = lambda_offset(i_layer + 1);
                    index_2 = lambda_offset(i_layer + 2);
                    for beta = 0 : index_1 - 1
                        P(M_index,beta + index_1, l_index) = sign(P(M_index,beta + index_2, l_index)) *...
                            sign(P(M_index,beta + index_1 + index_2, l_index)) * min(abs(P(M_index,beta + index_2, l_index)), abs(P(M_index,beta + index_1 + index_2, l_index)));
                    end
                end
            case N/2
                index_1 = lambda_offset(m);
                for beta = 0 : index_1 - 1
                    x_tmp = C(M_index,beta + index_1, 2 * l_index - 1);
                    P(M_index,beta + index_1, l_index) = (1 - 2 * x_tmp) * llr(M_index,beta + 1) + llr(M_index,beta + 1 + index_1);
                end
                for i_layer = m - 2 : -1 : 0
                    index_1 = lambda_offset(i_layer + 1);
                    index_2 = lambda_offset(i_layer + 2);
                    for beta = 0 : index_1 - 1
                        P(M_index,beta + index_1, l_index) = sign(P(M_index,beta + index_2, l_index)) *...
                            sign(P(M_index,beta + index_1 + index_2, l_index)) * min(abs(P(M_index,beta + index_2, l_index)), abs(P(M_index,beta + index_1 + index_2, l_index)));
                    end
                end
            otherwise
                index_1 = lambda_offset(layer + 1);
                index_2 = lambda_offset(layer + 2);
                for beta = 0 : index_1 - 1
                    P(M_index,beta + index_1, l_index) = (1 - 2 * C(M_index,beta + index_1, 2 * l_index - 1)) * P(M_index,beta + index_2, lazy_copy(M_index,layer + 2, l_index)) +...
                        P(M_index,beta + index_1 + index_2, lazy_copy(M_index,layer + 2, l_index));
                end
                for i_layer = layer - 1 : -1 : 0
                    index_1 = lambda_offset(i_layer + 1);
                    index_2 = lambda_offset(i_layer + 2);
                    for beta = 0 : index_1 - 1
                        P(M_index,beta + index_1, l_index) = sign(P(M_index,beta + index_2, l_index)) *...
                            sign(P(M_index,beta + index_1 + index_2, l_index)) * min(abs(P(M_index,beta + index_2, l_index)),...
                            abs(P(M_index,beta + index_1 + index_2, l_index)));
                    end
                end
        end
    end
    end
    
%%
    if frozen_bits(phi + 1) == 0%if now we decode an unfrozen bit
%% PM updating
        PM_pair = realmax * ones(M,2, L);

        for M_index=1:M
        for l_index = 1 : L
            if activepath(M_index,l_index) == 0
                continue;
            end
            if P(M_index,1, l_index) >= 0
                PM_pair(M_index,1, l_index) = PM(M_index,l_index);
                PM_pair(M_index,2, l_index) = PM(M_index,l_index) + P(M_index,1, l_index);
            else
                PM_pair(M_index,1, l_index) = PM(M_index,l_index) - P(M_index,1, l_index);
                PM_pair(M_index,2, l_index) = PM(M_index,l_index);
            end
        end
        end
       compare=zeros(M,2,L);
        for M_index=1:M
        middle = min(2 * sum(activepath(M_index,:)), L);
        [~,PM_inx] = sort(PM_pair(M_index,:));
        compare(M_index,PM_inx(1:middle))=1; 
        end
        
        
        kill_index = zeros(M,L, 1);%to record the index of the path that is killed
        kill_cnt = zeros(M,1);%the total number of killed path
       for M_index=1:M
         for i = 1 : L
            if (compare(M_index,1, i) == 0)&&(compare(M_index,2, i) == 0)%which indicates that this path should be killed
                activepath(M_index,i) = 0;
                kill_cnt(M_index) = kill_cnt(M_index) + 1;%push stack
                kill_index(M_index,kill_cnt(M_index),1) = i;
            end
         end
       end
%         for i=1:M
%             activepath(i,(compare(M,1, :) == 0)&(compare(M,2, :) == 0))=0;
%                 kill_cnt(i) = kill_cnt(i) + 1;%push stack
%                 kill_index(i,1:kill_cnt(i)) = find(activepath(i,:)==0);
%         end
        
       %% 路径判决
       for M_index=1:M
        for l_index = 1 : L
            if activepath(M_index,l_index) == 0
                continue;
            end
            path_state = compare(M_index,1, l_index) * 2 + compare(M_index,2, l_index);
            switch path_state%path_state can equal to 0, but in this case we do no operation.
                case 1
                    u(M_index,cnt_u, l_index) = 1;
                    C(M_index,1, 2 * l_index - 1 + phi_mod_2) = 1;
                    PM(M_index,l_index) = PM_pair(M_index,2, l_index);
                case 2
                    u(M_index,cnt_u, l_index) = 0;
                    C(M_index,1, 2 * l_index - 1 + phi_mod_2) = 0;
                    PM(M_index,l_index) = PM_pair(M_index,1, l_index);
                case 3
                    index = kill_index(M_index,kill_cnt(M_index),1);
                    kill_cnt(M_index) = kill_cnt(M_index) - 1;%pop stack
                    activepath(M_index,index) = 1;
                    %lazy copy
                    lazy_copy(M_index,:, index) = lazy_copy(M_index,:, l_index);
                    u(M_index,:, index) = u(M_index,:, l_index);
                    u(M_index,cnt_u, l_index) = 0;
                    u(M_index,cnt_u, index) = 1;
                    C(M_index,1, 2 * l_index - 1 + phi_mod_2) = 0;
                    C(M_index,1, 2 * index - 1 + phi_mod_2) = 1;
                    PM(M_index,l_index) = PM_pair(M_index,1, l_index);
                    PM(M_index,index) = PM_pair(M_index,2, l_index);
            end
        end
       end
        cnt_u = cnt_u + 1;
        %% 非冻结bit
    else%frozen bit operation
        for M_index=1:M
        for l_index = 1 : L
            if activepath(M_index,l_index) == 0
                continue;
            end
            if (P(M_index,1, l_index) < 0) &&(polar_info_esti_all(M_index,phi+1)==0)
                PM(M_index,l_index) = PM(M_index,l_index) - P(M_index,1, l_index);
            end
            if (P(M_index,1, l_index) > 0) &&(polar_info_esti_all(M_index,phi+1)==1)
                PM(M_index,l_index) = PM(M_index,l_index) + P(M_index,1, l_index);
            end
            if phi_mod_2 == 0
                C(M_index,1, 2 * l_index - 1) = polar_info_esti_all(M_index,phi+1);
            else
                C(M_index,1, 2 * l_index) = polar_info_esti_all(M_index,phi+1);
            end 
        end
        end
    end 
    
    for M_index=1:M
    for l_index = 1 : L%partial-sum return
        if activepath(M_index,l_index) == 0
            continue
        end
        if (phi_mod_2  == 1) && (phi ~= N - 1)
            layer = bit_layer_vec(phi + 1);
            for i_layer = 0 : layer - 1
                index_1 = lambda_offset(i_layer + 1);
                index_2 = lambda_offset(i_layer + 2);
                for beta = index_1 : 2 * index_1 - 1
                    C(M_index,beta + index_1, 2 * l_index) = mod(C(M_index,beta, 2 *  lazy_copy(M_index,i_layer + 1, l_index) - 1) + C(M_index,beta, 2 * l_index), 2);%Left Column lazy copy
                    C(M_index,beta + index_2, 2 * l_index) = C(M_index,beta, 2 * l_index);   
                end
            end
            index_1 = lambda_offset(layer + 1);
            index_2 = lambda_offset(layer + 2);
            for beta = index_1 : 2 * index_1 - 1
                C(M_index,beta + index_1, 2 * l_index - 1) = mod(C(M_index,beta, 2 * lazy_copy(M_index,layer + 1, l_index) - 1) + C(M_index,beta, 2 * l_index), 2);%Left Column lazy copy
                C(M_index,beta + index_2, 2 * l_index - 1) = C(M_index,beta, 2 * l_index);
            end 
        end
    end
    end
    %lazy copy
    if phi < N - 1
        for i_layer = 1 : llr_layer_vec(phi + 2) + 1
            
            for l_index = 1 : L
                lazy_copy(:,i_layer, l_index) = l_index;
            end
        end
    end
end
%path selection.
[~, path_ordered] = sort(PM,2);
polar_info_esti=zeros(M,K);
for M_index=1:M
polar_info_esti(M_index,:) = u(M_index,:, path_ordered(M_index,1));
end
j=1;
for i=1:N
if frozen_bits(i)==0
    polar_info_esti_all(:,i)=polar_info_esti(:,j);
    j=j+1;
end
end
end