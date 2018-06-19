clear all
clc
close all
tic
%参数设定
M=10;
DOA =[30 45 60]*pi/180;
SNR=10;
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
    w = (inv(R)*a)/(a'*inv(R)*a);
    y(i) = sqrt(abs(a'*R*a));
    y2(i) = sqrt(abs(w'*R*w));
end
doa = ESA(angle, y, P)
doa1 = ESA(angle, y2, P)
plot(angle,20*log(y/max(y)),angle,20*log(y2/max(y2)))
xlabel('\theta/degree');
ylabel('Sgnal Magnitude/dB');
legend('CBF算法','MVDR算法','Location','NorthWest')
toc