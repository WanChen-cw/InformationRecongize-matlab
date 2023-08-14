% 二进制删除信道的极化码构造(以信道容量为衡量标准）
%BEC信道巴氏参数计算  可靠性估计
function ZWi=BEC(N,epsilon)
% N=8;
% epsilon=0.5;  % 删除概率
n=log2(N);
ZWi=zeros(n+1,N);
ZWi(1,1)=epsilon;
for i=1:n
    k=2^(i-1);
    for j=1:k
        tmp=ZWi(i,j);
        ZWi(i+1,2*j)=tmp^2;
        ZWi(i+1,2*j-1)=2*tmp-tmp^2;  
    end 
end
scatter((1:N),ZWi(n+1,1:N),'.b');
axis([0 N 0 1]);
xlabel('Channel index');
ylabel('zw');
end