function[vHat,n] = decodeLLR_BP(llr, H, iteration,s)
%iteration 最大译码次数
%s 发送方校验子
%
%
[M,N] = size(H);

% 先验概率对数似然值
Lci = llr;

%E校验矩阵Lr
Lrji = zeros(M, N);
% M变量矩阵初始化
Lqij = H.*repmat(Lci, M, 1);

n=0;
while n <iteration
    %译码轮次
   fprintf('Iteration : %d\n', n);
   
   
   
    alphaij = sign(Lqij);  
    betaij = abs(Lqij);
   % 校验矩阵更新
   for i = 1:M
      % Find non-zeros in the column
      c1 = find(H(i, :));
      %符号信息
      if s(i)==1
          sg=-1;
      else
          sg=1;
      end
      for k = 1:length(c1)  
          c2=c1;
          c2(k)=[];
          Lrji(i, c1(k))=sg*prod(alphaij(i,c2))*min(betaij(i,c2));
%         prodofdate=prod(tanh(Lqij(i,c1)/2))/tanh(Lqij(i,c1(k))/2);
%         Lrji(i, c1(k))=sg*2*atanh(prodofdate);
      end % for k 
   end % for i
   
   % 变量矩阵更新
   
   LQ=zeros(1,N);
   for j = 1:N

      % Find non-zero in the row
      r1 = find(H(:, j));
      for k = 1:length(r1)         
         
         % Update L(qij) by summation of L(rij)\r1(k)
         Lqij(r1(k), j) = Lci(j) + sum(Lrji(r1, j)) - Lrji(r1(k), j);
      
      end % for k
      
      % 后验概率
      LQ(j) = Lci(j) + sum(Lrji(r1, j));
    end % for j 
      % 硬判决
      
        vHat=(LQ<0);
   %与Alice校验值比对
   if isequal(mod(vHat*H',2), s)
       break ;
    else
        n=n+1 ;
    end
   
end % for n

