N=2^3;
G=zeros(2^3,2^3);
for i=1:N
for j=1:N
    G(i,j)=GG(i,j,N);
end
end