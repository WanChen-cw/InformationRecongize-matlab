function [Ksa,Ksb]=provide(N,Eu)
    Ksa=rand(1,N) > 0.5;
    M=zeros(1,N);%alice筛选码
    M(1:floor(N*Eu))=1;
    K=randperm(N);
    Ksb=mod(Ksa+M(K),2);%此处ksb指bob筛选码
end