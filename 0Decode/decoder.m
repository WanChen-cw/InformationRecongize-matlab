    function polar_info_esti_all = decoder(llr, Sd)
    %输入为


    K=sum(Sd==-1);%信息比特数量
    frozen_bits=(Sd~=-1);%冻结bit位置
    llr = bitrevorder(llr);
    N = length(llr);
    polar_info_esti_all=Sd;


    %译码decode参数
    %分段向量，长度位n+1，元素值为1，2，4，8.。。   给存储中间变量结果的P C向量分段
    lambda_offset = 2.^(0 : log2(N));
    %llr计算实际执行层数。
    llr_layer_vec = get_llr_layer(N);
    bit_layer_vec = get_bit_layer(N);
    n = log2(N);
    P = zeros(N - 1, 1);%channel llr is not include in P.
    C = zeros(N - 1, 2);%C stores internal bit values
    polar_info_esti = zeros(K, 1);
    cnt_K = 1;
    for phi = 0 : N - 1
        switch phi
            case 0%for decoding u_1
                index_1 = lambda_offset(n);%index_1=N/2
                for beta = 0 : index_1 - 1%use llr vector
                    % f运算近似表达
                    P(beta + index_1) =  sign(llr(beta + 1)) * sign(llr(beta + 1 + index_1)) * min(abs(llr(beta + 1)), abs(llr(beta + 1 + index_1)));
                end
                for i_layer = n - 2 : -1 : 0 %use P vector
                    index_1 = lambda_offset(i_layer + 1);
                    index_2 = lambda_offset(i_layer + 2);
                    for beta = index_1 : index_2 - 1
                        P(beta) =  sign(P(beta + index_1)) * sign(P(beta + index_2)) * min(abs(P(beta + index_1)), abs(P(beta + index_2)));
                    end
                end
            case N/2%for deocding u_{N/2 + 1}
                index_1 = lambda_offset(n);
                for beta = 0 : index_1 - 1%use llr vector. g function.
                    P(beta + index_1) = (1 - 2 * C(beta + index_1, 1)) * llr(beta + 1) + llr(beta + 1 + index_1);
                end
                for i_layer = n - 2 : -1 : 0%use P vector. f function
                    index_1 = lambda_offset(i_layer + 1);
                    index_2 = lambda_offset(i_layer + 2);
                    for beta = index_1 : index_2 - 1
                        P(beta) =  sign(P(beta + index_1)) * sign(P(beta + index_2)) * min(abs(P(beta + index_1)), abs(P(beta + index_2)));
                    end
                end
            otherwise
                llr_layer = llr_layer_vec(phi + 1);
                index_1 = lambda_offset(llr_layer + 1);
                index_2 = lambda_offset(llr_layer + 2);
                for beta = index_1 : index_2 - 1%g function is first implemented.
                    P(beta) = (1 - 2 * C(beta, 1)) * P(beta + index_1) + P(beta + index_2);
                end
                for i_layer = llr_layer - 1 : -1 : 0%then f function is implemented.
                    index_1 = lambda_offset(i_layer + 1);
                    index_2 = lambda_offset(i_layer + 2);
                    for beta = index_1 : index_2 - 1
                        P(beta) =  sign(P(beta + index_1)) * sign(P(beta + index_2)) * min(abs(P(beta + index_1)), abs(P(beta + index_2)));
                    end
                end
        end
        phi_mod_2 = mod(phi, 2);
        if frozen_bits(phi + 1) == 1%frozen bit

            C(1, 1 + phi_mod_2) = polar_info_esti_all(phi+1);
        else%information bit
            C(1, 1 + phi_mod_2) = P(1) < 0;%store internal bit values
            polar_info_esti(cnt_K) = P(1) < 0;
            cnt_K = cnt_K + 1;
        end
        if phi_mod_2  == 1 && phi ~= N - 1
            bit_layer = bit_layer_vec(phi + 1);
            for i_layer = 0 : bit_layer - 1%give values to the 2nd column of C
                index_1 = lambda_offset(i_layer + 1);
                index_2 = lambda_offset(i_layer + 2);
                for beta = index_1 : index_2 - 1
                    C(beta + index_1, 2) = mod(C(beta, 1) + C(beta, 2), 2);
                    C(beta + index_2, 2) = C(beta, 2);
                end
            end
            index_1 = lambda_offset(bit_layer + 1);
            index_2 = lambda_offset(bit_layer + 2);
            for beta = index_1 : index_2 - 1%give values to the 1st column of C
                C(beta + index_1, 1) = mod(C(beta, 1) + C(beta, 2), 2);
                C(beta + index_2, 1) = C(beta, 2);
            end
        end
    end



    j=1;
    for i=1:N
    if frozen_bits(i)==0
        polar_info_esti_all(i)=polar_info_esti(j);
        j=j+1;
    end
    end

    end
