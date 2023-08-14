function [decode_date,cnt]=Bob(Ksb,N,frozen_bits,Pe,Eu,s1,s2,s3)
%译码函数

%% 输入polar译码分段
llr=zeros(1,length(Ksb));
llr(Ksb==0)=log2((1-Eu)/Eu);
llr(Ksb==1)=log2(Eu/(1-Eu));
number_polar=length(Ksb)/N;
llr=reshape(llr,N,number_polar)';
%% bit翻转
%llr=bitrevorder(llr')';
%% polar译码前置
[~,W]=sort(Pe, 'descend');
%初始化译码结果
decode_date=-1*ones(number_polar,N);
decode_date(:,W(1:length(frozen_bits(1,:))))=frozen_bits;
site_frozen_bits=(decode_date(1,:)~=-1);%冻结bit位置
number_frozen_bits=sum(site_frozen_bits);
%译码decode参数
%分段向量，长度位n+1，元素值为1，2，4，8.。。   给存储中间变量结果的P C向量分段
lambda_offset = 2.^(0 : log2(N));
%llr计算实际执行层数。
llr_layer_vec = get_llr_layer(N);
bit_layer_vec = get_bit_layer(N);
n = log2(N);
%中间存储
P = zeros(number_polar,N - 1);%channel llr is not include in P.
C = zeros(number_polar,N - 1, 2);%C stores internal bit values

%% ldpc 译码前置
H1 = IEEE80216e(number_polar, '1/2');
H1=logical(im2double(H1));
H2 = IEEE80216e(number_polar, '2/3A');
H2=logical(im2double(H2));
H3 = IEEE80216e(number_polar, '5/6');
H3=logical(im2double(H3));
number_ldpc=size(s1,1);%ldcp译码器个数
site_dopc=sort(W(number_frozen_bits+1:number_frozen_bits+number_ldpc));%ldpc译码位置
site_dopc(number_ldpc+1)=N+1;
K=1;
%% 联合译码
cnt=0;
for phi = 0 : N - 1
    switch phi
        case 0%for decoding u_1
            index_1 = lambda_offset(n);%index_1=N/2
            for beta = 0 : index_1 - 1%use llr vector
                % f运算近似表达
                P(:,beta + index_1) =  sign(llr(:,beta + 1)) .* sign(llr(:,beta + 1 + index_1)) .* min(abs(llr(:,beta + 1)), abs(llr(:,beta + 1 + index_1)));
            end
            for i_layer = n - 2 : -1 : 0 %use P vector
                index_1 = lambda_offset(i_layer + 1);
                index_2 = lambda_offset(i_layer + 2);
                for beta = index_1 : index_2 - 1
                    P(:,beta,1) =  sign(P(:,beta + index_1)) .* sign(P(:,beta + index_2)) .* min(abs(P(:,beta + index_1)), abs(P(:,beta + index_2)));
                end
            end
        case N/2%for deocding u_{N/2 + 1}
            index_1 = lambda_offset(n);
            for beta = 0 : index_1 - 1%use llr vector. g function.
                P(:,beta + index_1) = (1 - 2 .* C(:,beta + index_1, 1)) .* llr(:,beta + 1) + llr(:,beta + 1 + index_1);
            end
            for i_layer = n - 2 : -1 : 0%use P vector. f function
                index_1 = lambda_offset(i_layer + 1);
                index_2 = lambda_offset(i_layer + 2);
                for beta = index_1 : index_2 - 1
                    P(:,beta) =  sign(P(:,beta + index_1)) .* sign(P(:,beta + index_2)) .* min(abs(P(:,beta + index_1)), abs(P(:,beta + index_2)));
                end
            end
        otherwise
            llr_layer = llr_layer_vec(phi + 1);
            index_1 = lambda_offset(llr_layer + 1);
            index_2 = lambda_offset(llr_layer + 2);
            for beta = index_1 : index_2 - 1%g function is first implemented.
                P(:,beta,1) = (1 - 2 * C(:,beta, 1)) .* P(:,beta + index_1) + P(:,beta + index_2);
            end
            for i_layer = llr_layer - 1 : -1 : 0%then f function is implemented.
                index_1 = lambda_offset(i_layer + 1);
                index_2 = lambda_offset(i_layer + 2);
                for beta = index_1 : index_2 - 1
                    P(:,beta) =  sign(P(:,beta + index_1)) .* sign(P(:,beta + index_2)) .* min(abs(P(:,beta + index_1)), abs(P(:,beta + index_2)));
                end
            end
    end
    phi_mod_2 = mod(phi, 2);
    if site_frozen_bits(1,phi + 1) == 1%frozen bit 
        C(:,1, 1 + phi_mod_2) = decode_date(:,phi+1);
        if~isequal((P(:,1)<0), decode_date(:,phi+1))
            cnt=cnt+number_polar;
        end
    else%information bit
        if site_dopc(K)==(phi+1)
            if Pe(phi+1)<=0.07&&Pe(phi+1)>0.037
                    [vHat,item] = decodeLLR_BP(P(:,1)', H1, 40,s1(K,:));
                    if ~isequal(mod((P(:,1)'<0)*H1',2), s1(K,:))
                        cnt=cnt+length(s1(K,:));
                    end
            elseif Pe(phi+1)<=0.37&&Pe(phi+1)>0.013
                    [vHat,item] = decodeLLR_BP(P(:,1)', H2, 40,s2(K,:));
                    if ~isequal(mod((P(:,1)'<0)*H2',2), s2(K,:))
                        cnt=cnt+length(s2(K,:));
                    end  
            else
                    [vHat,item] = decodeLLR_BP(P(:,1)', H3, 40,s3(K,:));
                    if ~isequal(mod((P(:,1)'<0)*H3',2), s3(K,:))
                        cnt=cnt+length(s3(K,:));
                    end
            end

        
                C(:,1, 1 + phi_mod_2) =vHat;%store internal bit values
                decode_date(:,phi+1) = vHat;
                K=K+1;
        else
            C(:,1, 1 + phi_mod_2) =P(:,1)<0;%store internal bit values
            decode_date(:,phi+1) = P(:,1)<0;
        end
    end
    if phi_mod_2  == 1 && phi ~= N - 1
        bit_layer = bit_layer_vec(phi + 1);
        for i_layer = 0 : bit_layer - 1%give values to the 2nd column of C
            index_1 = lambda_offset(i_layer + 1);
            index_2 = lambda_offset(i_layer + 2);
            for beta = index_1 : index_2 - 1
                C(:,beta + index_1, 2) = mod(C(:,beta, 1) + C(:,beta, 2), 2);
                C(:,beta + index_2, 2) = C(:,beta, 2);
            end
        end
        index_1 = lambda_offset(bit_layer + 1);
        index_2 = lambda_offset(bit_layer + 2);
        for beta = index_1 : index_2 - 1%give values to the 1st column of C
            C(:,beta + index_1, 1) = mod(C(:,beta, 1) + C(:,beta, 2), 2);
            C(:,beta + index_2, 1) = C(:,beta, 2);
        end
    end
end




end