N=2^14;
Eu=0.02;
load('16384-0.02-PeD.mat')
C=-Eu*log2(Eu)-(1-Eu)*log2(1-Eu);
min=ceil(N*C);
number_frozen_bits=min:5:N;
[PD,W]=sort(PeD, 'descend');
fer_est_D=zeros(1,length(number_frozen_bits));

%%
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
semilogy(number_frozen_bits,fer_est_D);
title('码长 qber')
xlabel('冻结bit数量') 
ylabel('估计误帧率')
legend('降级信道失败率估计值')
name1='%d_%G误差估计与.fig';
name=sprintf(name1,N,Eu);
% saveas( gca, name)