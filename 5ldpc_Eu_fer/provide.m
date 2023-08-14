function [Ksa,Ksb]=provide(N,Eu)
    Ksa=rand(1,N) > 0.5;
    M=zeros(1,N);
    M(1:floor(N*Eu))=1;
    K=randperm(N);
    Ksb=mod(Ksa+M(K),2);
end