semilogy(f,fer_ldpc_polar,f_ir,fer_IR_polar,f_ldpc,fer_ldpc);
title(' ')
xlabel('f译码效率') 
ylabel('误帧率')
legend('级联误帧率','IR误帧率','ldpc')
% xlim([f(21) f(3)])
% name1='%d_%G对比结果.fig';
% % name=sprintf(name1,N,Eu);
saveas( gca, 'f对比图-2')