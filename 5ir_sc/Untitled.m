%传统polar应用模式性能分析——IR
N=2^16;
Eu=0.02;
load('65536_0.02.mat')
C=-Eu*log2(Eu)-(1-Eu)*log2(1-Eu);
min=ceil(N*C);
number_frozen_bits=10000:500:12000;
[PD,W]=sort(PeD, 'descend');

fer_est_D=zeros(1,length(number_frozen_bits));

 for loop=1:length(number_frozen_bits)

ber_D=PeD(sort(W(number_frozen_bits(loop)+1:N)));
frist_D=zeros(1,length(ber_D));
frist_D(1)=ber_D(1);
for i=2:length(ber_D)
frist_D(i)=prod(1-ber_D(1:i-1))*ber_D(i);
end
fer_est_D(loop)=sum(frist_D);
 end
 
 
 %% 画图
semilogy(number_frozen_bits(1:length(fer_est_D)),fer_est_D);
title('冻结bit数量对的影响')
xlabel('冻结bit数量') 
ylabel('译码失败概率')
legend('降级信道失败率估计值')
