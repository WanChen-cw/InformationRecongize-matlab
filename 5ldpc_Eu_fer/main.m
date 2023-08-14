%main
%N 码长    Eu 误码率    H 校验矩阵
N=1056;

H = IEEE80216e(1056, '5/6');
H=logical(im2double(H));
Eu=0.01:0.003:0.037;
mea_ldpc_R4=zeros(1,length(Eu));
for j=1:length(Eu)
mea_ldpc=0;
%%循环测试

for n=1:100


%数据生成
[Ksa,Ksb]=provide(N,Eu(j));
%先验概率对数似然值
llr=zeros(1,length(Ksb));
llr(Ksb==0)=log2((1-Eu(j))/Eu(j));
llr(Ksb==1)=log2(Eu(j)/(1-Eu(j)));
%检验子
s=rem(Ksa*H',2);
%译码

    tic;
    
    
[u_decode,n3]= decodeLLR_BP_sparse(llr,H, 50,s);

time = toc;
fprintf('代码运行时间为：%f 秒\n', time);


%译码成功计数
if isequal(u_decode,Ksa)
mea_ldpc=mea_ldpc+1;
end
end

mea_ldpc_R4(j)=mea_ldpc/10000;
end
plot(Eu,mea_ldpc_R4);
title('码长1056，码率1/2')
xlabel('Eu') 
ylabel('ldpc译码成功率')
% name1='%d_%G误差估计与实测.fig';
% name=sprintf(name1,N,Eu);
% saveas( gca, 'ldpc性能测试2/3.fig')