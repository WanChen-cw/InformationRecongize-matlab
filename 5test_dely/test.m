%delay test
N=2^18;
Eu=0.02;
%load('262144_0.02.mat')
%number_frozen_bits=50000;
number_frozen_bits=360;
[PD,W]=sort(Pe, 'descend');
%%循环测试
M=1;
[Ksa,Ksb]=provide(N*M,Eu);
M=length(Ksa)/N;
date_polar=reshape(Ksa,N,M)';
u=encode(Ksa,N);
frozen_bits=u(:,W(1:number_frozen_bits));
tic
u_est=decode_polar_m(Ksb,N,frozen_bits,Pe,Eu);
toc

disp(['运行时间: ',num2str(toc)]);
number=sum(sum(u~=u_est,2)~=0)/M;

