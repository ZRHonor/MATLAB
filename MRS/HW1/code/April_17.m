clear all
clc
tic
%参数设定
M=10;
DOA =[5 45 65]*pi/180;
SNR=10;
d = 0.5;
N = 400;
QAM = 16;
% X = StatSigGenerate(M, N, DOA, SNR, 'Independent', d);
P = length(DOA);
A=exp(-j*2*pi*0.5*[0:M-1]'*sin(DOA));
%信源模型建立
for k=1:P
    symbol = randi([0, QAM-1], 1, N);
    S(k,:) = qammod(symbol, QAM);
%     S(k,:)=sqrt(10.^(SNR/10))*sind(1:N);
end
temp = abs(S).^2;
signalpower = mean(temp(:));
X = awgn(A*S, SNR, signalpower);

%协方差矩阵特征值分解得到噪声子空间
R=X*X'/N;
angle = -90:0.01:90;
for i =1:length(angle)
	a = exp(-j*2*pi*0.5*[0:M-1]'*sin(pi*angle(i)/180));
	y(i) = sqrt(abs(a'*R*a));
end
plot(angle,y)
doa = ESA(angle, y, length(DOA))
toc