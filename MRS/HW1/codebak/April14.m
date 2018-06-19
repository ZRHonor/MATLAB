clear all;
clc;
close all;

M = 10;                         % 阵元数目
DOA = [5 45 60].*pi/180;		% 来波方向
d = 0.5;						% 阵元距离
N = 400;						% 采样点数
QAM = 16;						% 16QAM调制
SNR = 10;

% 阵列流型A
P = length(DOA);
for i = 1:P
	 A(:,i)=exp(-j*2*pi*[0:M-1]'*sin(DOA(i))*0.5);
end

% 信源模型
% TODO QAM调制？？？？？？？？？？？？？？？？？
for j = 1:P
	s = sind(1:N);
	S(j,:) = s;
end

% 接收信号模型建立
 X=awgn(A*S, SNR);

% 协方差矩阵
R=X*X'/N;

angle = -90:1:90;
for k = 1:length(angle)
    a = exp(-j*2*pi*0.5*[0:M-1]'*sin(angle(k)));
    y(k) = sqrt(abs(a'*R*a));
end
plot(angle, 10*log(y/max(y)))
axis([-90 90 -80 10]);
xlabel('theta/degree');
ylabel('归一化空间谱/dB');