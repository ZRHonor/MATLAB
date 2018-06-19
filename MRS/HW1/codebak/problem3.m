clear all
clc
close all
tic
%参数设定
M = 10;
DOA = [30 45 60]*pi/180;
SNR = 10;
d = 0.5;
N = 4000;
QAM = 16;
P = length(DOA);
A=exp(-j*2*pi*0.5*[0:M-1].'*sin(DOA));
%信源模型建立
for k=1:P
    symbol = randi([0, QAM-1], 1, N);
%     S(k,:) = qammod(symbol, QAM);
    S(k,:) = qammod(symbol, QAM);
end
% X = A*S;
scatterplot(S(2,:));
title('(a)原信号')
X = awgn(A*S,SNR,'measured');
% for i = 1:M
%     X(i,:) = awgn(X(i,:),SNR,'measured');
% end
% 阵元1接收到的信号
re_sig1 = X(1,:);
scatterplot(re_sig1);
title('(b)阵元1接收到的信号')

% 空间匹配滤波
theta = 45*pi/180;
a = exp(-j*2*pi*d*[0:M-1].'*sin(theta));

scatterplot(a'*X/M);
title('(a)CBF算法得到的信号')
% hold on
% MVDR
R = X*X'/N;
w = (inv(R)*a)/(a'*inv(R)*a);

scatterplot(w'*X);
title('(b)MVDR算法得到的信号')

% R1 = cov(X');
% w1 = (inv(R1)*a)/(a'*inv(R1)*a);
% scatterplot(w1'*X)
% title('MVDR2')

% ans1 = a'*X/M;
% ans2 = w'*X;
% 
% scatter(real(ans1),imag(ans1),1,'k')
% hold on
% scatter(real(ans2),imag(ans2),1,'r')