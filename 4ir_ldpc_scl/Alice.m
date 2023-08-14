function [frozen_bits,s,u]=Alice(Ksa,N,PeD,H,number_frozen_bits,number_ldpc)
% Ksa-输入数据   N-极化码码长   
% Eu-qber信道误码率    G-极化码生成矩阵
% Pe-极化码信道误码率   H-ldpc校验矩阵  numbe_ldpc-ldpc个数
%发送端协商后密钥 u

%% 输入分段并行极化编码
M=length(Ksa)/N;
GN=G(N);%极化码生成矩阵
[~,W]=sort(PeD, 'descend');
%C=-Eu*log2(Eu)-(1-Eu)*log2(1-Eu);
%number_frozen_bits=N-floor(N*C);
date_polar=reshape(Ksa,N,M)';
u=date_polar*GN;
u=mod(u,2);
frozen_bits=u(:,W(1:number_frozen_bits));
%% ldpc校验子生成
    date_ldpc=(u(:,sort(W(number_frozen_bits+1:number_frozen_bits+number_ldpc))))';
    s=rem(date_ldpc*H',2);
end