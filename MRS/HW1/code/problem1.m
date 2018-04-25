clear all
clc
close all
tic
%参数设定
M=10;
DOA =[5 45 65]*pi/180;
SNR=3;
d = 0.5;
N = 400;
QAM = 16;

% 建立信号模型
P = length(DOA);
A=exp(-j*2*pi*0.5*[0:M-1]'*sin(DOA));
for k=1:P
    symbol = randi([0, QAM-1], 1, N);
    S(k,:) = qammod(symbol, QAM);
end

X = awgn(A*S,SNR,'measured');

R=X*X'/N;
angle = -90:0.01:90;
for i =1:length(angle)
    a = exp(-j*2*pi*d*[0:M-1]'*sin(pi*angle(i)/180));
    y(i) = sqrt(abs(a'*R*a));
end
doa = ESA(angle, y, P)
plot(angle,10*log(y/max(y)))
xlabel('theta/degree');
ylabel('归一化空间谱/dB');
title('空间谱');
toc