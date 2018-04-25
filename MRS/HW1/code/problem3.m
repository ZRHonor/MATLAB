clear all
clc
close all
tic
%参数设定
M = 10;
DOA = [5 45 65]*pi/180;
SNR = 10;
d = 0.5;
N = 400;
QAM = 16;
P = length(DOA);
A=exp(-j*2*pi*0.5*[0:M-1]'*sin(DOA));
%信源模型建立
for k=1:P
    symbol = randi([0, QAM-1], 1, N);
%     S(k,:) = qammod(symbol, QAM);
    S(k,:) = qammod(symbol, QAM);
end
X = awgn(A*S,SNR,'measured');
% 阵元1接收到的信号
re_sig1 = X(1,:);
scatterplot(re_sig1);

% 空间匹配滤波
theta = 45*pi/180;
a = exp(-j*2*pi*d*[0:M-1]'*sin(theta));
scatterplot(a'*X/M);
title('空间匹配滤波')
% MVDR
R = X*X'/N;
w = inv(R)*a/(a'*inv(R)*a);
scatterplot(w'*X);
title('MVDR')

