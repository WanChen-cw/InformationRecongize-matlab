%生成矩阵
function GN = G(N) 
n=log2(N);
G=1;
for i=1:n
G=kron(G,[1,0;1,1]);
end
%GN=bitrevorder(G);%执行bit反转
 GN=G;
end