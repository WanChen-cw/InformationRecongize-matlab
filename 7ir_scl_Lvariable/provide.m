function [Ksa,Ksb]=provide(N,Eu)
Ksa=rand(1,N) > 0.5;
    M=zeros(1,N);%alice筛选码
    M(1:floor(N*Eu))=1;
    K=randperm(N);
    Ksb=mod(Ksa+M(K),2);%此处ksb指bob筛选码
for i=1:N
if Ksb(i)==0
    Ksb(i)=log((1-Eu)/Eu);%此处KSb指对数似然比
else 
    Ksb(i)=log(Eu/(1-Eu));
end
end