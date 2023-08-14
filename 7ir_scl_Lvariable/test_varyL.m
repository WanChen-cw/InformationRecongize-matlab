%optimized_frozen_vectors
N=2^16;
Eu=0.02;
t=10000;%误帧率测试轮次
Lmax=1;
load('65536_0.02.mat')
[~,W]=sort(Pe, 'descend');
n_fro=10000:500:12000;
result=zeros(length(n_fro),t,3);

for loop=1:length(n_fro)
    for i=1:t
    L=1;
    [Ksa,Ksb]=provide(N,Eu);
    u=encode(Ksa,N);
    Sd=-1*ones(1,N);
    S=Index_select(u,W(1:n_fro(loop)));
    Sd=Insert(Sd,S,W(1:n_fro(loop)));
    uu=SCL_decoder(Ksb,L,Sd);
%     
%     Sd16=-1*ones(1,N);
%     S16=Index_select(u,W(1:n_fro(loop)));
%     Sd16=Insert(Sd16,S16,W(1:n_fro(loop)));
%     uu16=SCL_decoder(Ksb,16,Sd16);
    while((L<Lmax)&&(~isequal(uu,u)))
        L=2*L;
        uu=SCL_decoder(Ksb,L,Sd);
    end
    result(loop,i,1)=isequal(uu,u);
    result(loop,i,2)=L;
%     result(loop,i,3)=isequal(uu16,u);
    end
end