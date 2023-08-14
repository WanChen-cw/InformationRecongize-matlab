%数据处理
clc
load('matlab.mat')
%1 list 1-32变化译码01表示失败成功
%2 list 表示译码最大list数量
%3 固定list为16时 0 1译码成功失败
fer1=1-sum(result(:,:,1),2)/length(result(1,:,1));
fer2=1-sum(result(:,:,3),2)/length(result(1,:,3));
fer3=1-sum(result(:,:,4),2)/length(result(1,:,3));
loop=sum(log2(result(:,:,2))+1,2)/length(result(1,:,3));
n_fro=10000:500:12000;
R=1-n_fro/(2^16);
f=(1-R)/0.1414;
figure;
% burg法
subplot(121);
plot(f,fer1,f,fer3);
xlabel('f因子')
xticks(f)
legend('list自适应,lmax=32','sc')
title('fer')
grid on;
subplot(122);
plot(f,loop);
title('平均译码轮次')
xticks(f)
xlabel('f因子')
grid on;

