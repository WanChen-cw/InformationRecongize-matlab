%输入损耗db=10*log10(1/q)

%发射信息长度L0，depletion损耗


% L0=100000000;
% depletion=30;
% q=1/(10^(depletion/10));
% [L,n]=optimze(q);
% bit_screen_length=L0*L;

depletion=1:80;
q=1./(10.^(depletion/10));
for i=1:80
L(i,:)=optimze(q(i));
end

plot(depletion,L(:,1));