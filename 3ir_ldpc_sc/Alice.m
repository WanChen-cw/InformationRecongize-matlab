function [frozen_bits,s1,s2,s3,u]=Alice(Ksa,N,PeD,number_frozen_bits,number_ldpc)
% Ksa-输入数据   N-极化码码长   
% Eu-qber信道误码率    G-极化码生成矩阵
% Pe-极化码信道误码率   H-ldpc校验矩阵  numbe_ldpc-ldpc个数
%发送端协商后密钥 u

%% 输入分段并行极化编码
M=length(Ksa)/N;
date_polar=reshape(Ksa,N,M)';

[~,W]=sort(PeD, 'descend');
%C=-Eu*log2(Eu)-(1-Eu)*log2(1-Eu);
%number_frozen_bits=N-floor(N*C);
u=zeros(M,N);
% 按位生成
% % G2=G(2^2);%极化码生成矩阵
% % load('G16384.mat', 'G14')
% % for i=1:4
% %     for j=1:4
% % u(1:M,1+(i-1)*2^14:(i)*2^14)=u(1:M,1+(i-1)*2^14:(i)*2^14)+date_polar(:,1+(j-1)*2^14:(j)*2^14)*kron(G2(j,i),G14);
% %     end
% % end
% % clear G14

GN=G(N);
u=date_polar*GN;
u=mod(u,2);
frozen_bits=u(:,W(1:number_frozen_bits));
%% ldpc校验子生成

H1 = IEEE80216e(M, '1/2');
H1=logical(im2double(H1));
H2 = IEEE80216e(M, '2/3A');
H2=logical(im2double(H2));
H3 = IEEE80216e(M, '5/6');
H3=logical(im2double(H3));
%     site_ldpc=sort(W(number_frozen_bits+1:number_frozen_bits+number_ldpc));
    date_ldpc=(u(:,sort(W(number_frozen_bits+1:number_frozen_bits+number_ldpc))))';
    % 实际传输每个ldpc只需要一个s
    s1=rem(date_ldpc*H1',2);
    s2=rem(date_ldpc*H2',2);
    s3=rem(date_ldpc*H3',2);
end