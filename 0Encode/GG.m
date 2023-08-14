function GN=GG(i,j,N)
%按位生成矩阵
n=log2(N);
x=logical(dec2bin(i-1,n)=='1');
y=logical(dec2bin (j-1,n)=='1');
GN=prod(xor(1,xor(y,(x&y))));
end