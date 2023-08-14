G0=G(2);
G3=G(1024);
G4=G(512);
G1=kron(G4,G0(1,:));
G2=kron(G4,G0(2,:));
daan=isequal(G3(1:512,:),G1);





% %编码运行时间对比
% u=rand(1,8)>0.5;
% N=length(u);
% tic
% G10= G(N);
% x1=mod(u*G10,2);
% toc
% disp(['运行时间: ',num2str(toc)]);
% 
% tic
% x2=zeros(1,N);
% for i=1:N
%     for j=1:N
% x2(i)=x2(i)+u(j)*GG(j,i,N);
%     end
% end
% x2=mod(x2,2);
% % G11= G1(N);
% % x2=mod(u*G11,2);
% 
% toc
% disp(['运行时间: ',num2str(toc)]);
% daan=isequal(x1,x2);



% %生成矩阵运行时间对比
% tic
% G10= G1(16);
% toc
% disp(['运行时间: ',num2str(toc)]);
% tic
% G11=G2(16);
% toc
% disp(['运行时间: ',num2str(toc)]);
% daan=isequal(G10,G11);




% %不同信道估计方法对比
% addpath('D:\studing\yanjiusheng\后处理\post-processing-m\2IR\AIR/')
% N=2^11;
% GN=G(N);
% Eu=0.02;
% %  Pe=BSC(N,Eu); 
% [~,W1]=sort(Pe, 'descend');
% 
% for j=1:1
% Ecnt1=0;
%    R=0.15;
% for i=1:10000
% [Ksa,Ksb]=provide(N,Eu);
% u=mod(Ksa*GN,2);
% Sd1=-1*ones(1,N);
% S1=Index_select(u,W1(1:ceil(N*R)));
% Sd1=Insert(Sd1,S1,W1(1:ceil(N*R)));
% u1=SCL_decoder(Ksb,16,Sd1);
% if ~isequal(u1,u)
%  Ecnt1=Ecnt1+1;
% end
% end
% p1(j)=Ecnt1/10000;
% end