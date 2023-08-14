function GN=G1(N)
%按位生成矩阵
n=log2(N);
GN=ones(N,N);
x=zeros(1,n);
y=zeros(1,n);
for i=1:N
    x=dec2bin(i-1,n);
    for j=1:N
        y=dec2bin (j-1,n);
        for m=1:n
        GN(i,j)=GN(i,j)*mod(1+str2num(y(m))+str2num(x(n+1-m))*str2num(y(m)),2);
        end
    end
end
% GN=bitrevorder(GN);
end