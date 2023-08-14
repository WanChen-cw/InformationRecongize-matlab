clear
n = 10;
N = 2^n;
miu = 64;
v = miu/2;
W = [1-0.02, 0.02;0.02,1-0.02];
IW = get_BMS_capacity(W);
Pe = bit_channel_degrading_procedure(W, 1 : N, miu);