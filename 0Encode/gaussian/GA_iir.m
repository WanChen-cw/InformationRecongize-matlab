%高斯信道对数似然值估计
%返回对数似然比
function iir=GA_iir(sigma,N)
% snr = 1;%信噪比
% R = 0.5;%码率
% N=512;%长度
iir = zeros(1, N);
iir(1) = 2/sigma^2;
% for i = 1:log2(N)
%     v=iir;
%     j = 2^(i - 1);
%     for k = 1:j
%         tmp = v(k);
%         iir(2*k-1) = phi_inverse(1 - (1 - phi(tmp))^2);
%         iir(2*k) = 2 * tmp;
%     end
% end
for i = 1:log2(N)
    j = 2^(i - 1);
    for k = 1:j
        tmp = iir(k);
        iir(k) = phi_inverse(1 - (1 - phi(tmp))^2);
        iir(k+j) = 2 * tmp;
    end
end
iir = bitrevorder(iir);
% scatter((1:N),iir(1:N),'.b');
% axis([0 1.1*N 0 4*N]);
% xlabel('Channel index');
% ylabel('E(LLRi)');
% toc
end