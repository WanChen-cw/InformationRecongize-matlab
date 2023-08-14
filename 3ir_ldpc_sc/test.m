%test ldpc——polar

%% 信道参数
N=2^10;%极化码码长
M=;%ldpc码长
Eu=0.02;%信道误码率（bsc转移概率）
 %load('65536_0.02.mat')
load('1024_0.02PeD.mat')
 %load('2048_0.02.mat')
%极化码信道可靠性估计

H1 = IEEE80216e(M, '1/2');
H1=logical(im2double(H1));
H2 = IEEE80216e(M, '2/3A');
H2=logical(im2double(H2));
H3 = IEEE80216e(M, '5/6');
H3=logical(im2double(H3));
% H1 = IEEE80216e(1056, '2/3A');0.37
% H1=logical(im2double(H1));
%1/2 0.7 2/3 0.37 5/6 0.13
%% bit分类
%ldpc性能参数 瀑布区门限
%冻结bit数量 根据极化信道转移概率 使转移概率适合与ldpc瀑布区门限
number_frozen_bits= sum(PeD>0.07);
%可靠bit数量门限_误帧率fer<10^-6   R=440/1024 误帧率由极化信道错误率估计 此值未考虑ldpc译码失败
%number_info=[300 370];%冻结bit下限
number_info=540;
% number_info=[14010 14330];
%ldcp_bit数量 ldpc码率固定时 可以译码的瀑布去门限
number_ldpc=number_info-number_frozen_bits;
%
f=(number_frozen_bits+number_ldpc*0.5)/(N*0.1414);
fer_mea=zeros(1,length(number_ldpc));

for i=1:length(number_ldpc)
number=0;
for loop=1:1000
[Ksa,Ksb]=provide(N*M,Eu);
[frozen_bits,s1,s2,s3,u]=Alice(Ksa,N,PeD,number_frozen_bits,number_ldpc(i));
tic
[decode_date,cnt]=Bob(Ksb,N,frozen_bits,PeD,Eu,s1,s2,s3);

toc
%disp(['运行时间: ',num2str(toc)]);
number=number+sum(sum(u~=decode_date,2)~=0);
end
fer_mea(i)=number/1000/M;
end