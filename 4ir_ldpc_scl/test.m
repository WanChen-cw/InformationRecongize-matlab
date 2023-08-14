%test ldpc——polar

%% 信道参数
N=1024;%极化码码长
M=1056;%ldpc码长
Eu=0.02;%信道误码率（bsc转移概率）
%极化码信道可靠性估计
load('1024_0.02PeD.mat')
H = IEEE80216e(1056, '1/2');
H=logical(im2double(H));

%% bit分类
%ldpc性能参数 瀑布区门限
%冻结bit数量 根据极化信道转移概率 使转移概率适合与ldpc瀑布区门限
C=-Eu*log2(Eu)-(1-Eu)*log2(1-Eu);
number_frozen_bits= sum(PeD>0.067);
%可靠bit数量门限_误帧率fer<10^-6   R=440/1024 误帧率由极化信道错误率估计 此值未考虑ldpc译码失败
number_info=200:10:400;%冻结bit下限
%ldcp_bit数量 ldpc码率固定时 可以译码的瀑布去门限
number_ldpc=number_info-number_frozen_bits;
%
f=(number_frozen_bits+number_ldpc*0.5)/(N*0.1414);
fer_mea=zeros(1,length(number_ldpc));
number=zeros(length(number_ldpc),10000);
for i=3:length(number_ldpc)

for loop=1:1000
[Ksa,Ksb]=provide(N*M,Eu);
[frozen_bits,s,u]=Alice(Ksa,N,PeD,H,number_frozen_bits,number_ldpc(i));
[decode_date]=Bob2(Ksb,N,16,frozen_bits,PeD,Eu,s,H);
number(i,loop)=sum(sum(u~=decode_date,2)~=0)/M;
end
end