function [L1,n]=optimze(q)
%输入响应率
%输出MZRLFL压缩效率、输出最优化n值

z=floor(-log2(-log(1-q)));
L1=L(z,q);
while(1)
z=z+1;
L2=L(z,q);
if L1<=L2
break;
else 
    L1=L2;
end
end
n=2^(z-1);
end

function L=L(z,q)
L=q*z/(1-(1-q)^(2^z-1));
end


% function [L1,n]=optimze(q,Nmax)
% %输入响应率
% %输出MZRLFL压缩效率、输出最优化n值
% %最优化n值受系统资源限制n<Nmax
% 
% z=floor(-log2(-ln(1-q)));
% L1=L(z,q);
% while(1)
% z=z+1;
% L2=L(z,q);
% if L1<=L2
% break;
% else 
%     L1=L2;
% end
% end
% n=2^(z-1);
% zmax=floor(log2(Nmax));
% if L(2^zmax)>L(Nmax)
%     n=Nmax;
%     L1=L(n,q);
% else
%     n=2^zmax;
%     L1=L(n,q);
% end
% 
% end