function x=encode(u,N)
v=1:N;
x=zeros(1,N);
for i=1:log2(N)
v1=reshape(v,2^i,[]);
v2=reshape(v1(1:2^(i-1),:),1,[]);
x(v2)=xor(u(v2),u(v2+2^(i-1)));
x(v2+2^(i-1))=u(v2+2^(i-1));
u=x;
end
end