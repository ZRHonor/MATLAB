clc; clear all; close all;
disp('你好nihao')
% 波速
% c = physconst('LightSpeed');
% % 频率
% fc = 300e6;
% % 波长
% lambda = c/fc;

% 线阵的阵元数
M = 10;
% 采样点数
N = 400;
% 阵元间距
d =  0.5;
% 信噪比
SNR = 10;
% 16QAM调制
QAM = 16;
% 信源方向
phi_s = [5, 45, 60];

P = length(phi_s);

% 阵列流型
for i=1:P
	% a 	= exp(-j*w*t)
	% 	= exp(-j*2*pi*f*(d*sin(theta)/(f*lambda)))
	% 	= exp(-j*2*pi*d*sin(theta)/lambda)
        A(:,i)=exp(-j*2*pi*[0:M-1]'*sin(phi_s(i))*0.5);
  end

% 信源模型
for k=1:P
	% 16QAM调制
	% 基带信号，不考虑载波
    	QAM = 16;
    	x = randi([0 QAM-1],1,N);
    	s = qammod(x, QAM);
    	S(k,:) = s;
end

% 阵元接收到的信号
X = awgn(A*S, SNR);

R=X*X'/N


angle = -90:1:90;
for i = 1:length(angle)
	a = exp(-j*2*pi*0.5*[0:M-1]'*sin(angle(i)));
	y(i) = sqrt(abs(a'*R*a))/3;
end
plot(angle,y)


% for i = 1:P
% 	u_s  = (d/lambda)*sin(phi_s(1)*pi/180);
% 	c_mf = exp(-i*2*pi*u_s*(0:M-1).')/sqrt(M);
% 	s = zeros(M,N);
% 	s(:,100) = 10^(SNR/20)*exp(-i*2*pi*u_s*(0:M-1).')/sqrt(M);
% 	w = (randn(M,N)+i*randn(M,N))/sqrt(2);
% 	x = s + w;
% end

% % Target Normalized Spatial Frequency.
% u_s  = (d/lambda)*sin(phi_s*pi/180);

% % 建立信号源模型
% P = length(phi_s)
% s = zeros( M, N, P);
% s(:, 200, :) = 10^(SNR/20)*exp(-i*2*pi*u_s'*((0:M-1)))/sqrt(M);

% % 高斯噪声
% w = (randn(M, P, N)+i*randn(M, P, N))/sqrt(2);

% c_mf = exp(-i*2*pi*u_s'*((0:M-1)))/sqrt(M);


