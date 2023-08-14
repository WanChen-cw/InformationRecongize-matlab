%传统polar应用模式性能分析——IR
N=2^16;
Eu=0.02;
C=-Eu*log2(Eu)-(1-Eu)*log2(1-Eu);
min=ceil(N*C);
number_frozen_bits=10000:500:12000;
[PD,W]=sort(PeD, 'descend');
[PU,W1]=sort(PeU, 'descend');

fer_est_D=zeros(1,length(number_frozen_bits));
fer_est_U=zeros(1,length(number_frozen_bits));
fer_mea=zeros(1,length(number_frozen_bits));
sigma2=zeros(1,length(number_frozen_bits));
number=zeros(length(number_frozen_bits),10000);
 for loop=1:length(number_frozen_bits)
ber_D=PeD(sort(W(number_frozen_bits(loop)+1:N)));
frist_D=zeros(1,length(ber_D));
frist_D(1)=ber_D(1);
for i=2:length(ber_D)
frist_D(i)=prod(1-ber_D(1:i-1))*ber_D(i);
end
fer_est_D(loop)=sum(frist_D);

ber_U=PeU(sort(W1(number_frozen_bits(loop)+1:N)));
frist_U=zeros(1,length(ber_U));
frist_U(1)=ber_U(1);
for i=2:length(ber_U)
frist_U(i)=prod(1-ber_U(1:i-1))*ber_U(i);
end
fer_est_U(loop)=sum(frist_U);  
%       fer_est_D_2(loop)=1-prod(1-PD(number_frozen_bits(loop)+1:N));
%       fer_est_U_2(loop)=1-prod(1-PU(number_frozen_bits(loop)+1:N));
 for i=1:1000  
%输入数据生成

M=10;
[Ksa,Ksb]=provide(N*M,Eu);

M=length(Ksa)/N;
date_polar=reshape(Ksa,N,M)';
u=date_polar*GN;
u=mod(u,2);
frozen_bits=u(:,W(1:number_frozen_bits(loop)));
tic
u_est=decode_polar_m(Ksb,N,frozen_bits,PeD,Eu);
toc

number(loop,i)=sum(sum(u~=u_est,2)~=0)/M;
% for j=1:M
%     number=number+isequal(u_est(j,:),u(j,:));
% end

end
fer_mea(loop)=mean(number(loop,:));
sigma2(loop)=var(number(loop,:));
if fer_mea(loop)==0
    break
end
 end
 
 
 %% 画图
semilogy(number_frozen_bits(1:length(fer_est_D)),fer_est_D,number_frozen_bits(1:length(fer_est_D)),fer_est_U,number_frozen_bits(1:length(fer_est_D)),fer_mea,number_frozen_bits(1:length(fer_est_D)),sigma2);
title('冻结bit数量对的影响')
xlabel('冻结bit数量') 
ylabel('译码失败概率')
legend('降级信道失败率估计值','升级信道失败率估计值','失败率测量值','方差')
name1='%d_%G误差估计与实测.fig';
name=sprintf(name1,N,Eu);
saveas( gca, name)