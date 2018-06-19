% 距离向处理

clear all; clc; close all;
%% 初始化
c=3e8;                                      % 光速 (m/s)
lambda = 1.5;                               % 波长 (m)
fc = c/lambda;                              % 载频 (Hz)
B = 1e6;                                    % 调频带宽 (Hz)
Tp = 400e-6;                                % 发射脉宽 (s)
PRT = 4000e-6;                              % 脉冲重复周期 (s)
M = 16;                                     % 阵元个数
d = 0.8;                                    % 阵元间距 (m)
PN = 8;                                     % pulse number

%% 参数设置
Fs = 2*B;
Ts = 1/Fs;
Rmin = 50e3; Rmax = 250e3;
Rm = Rmax - Rmin;

Nrn = ceil(2*Fs*Rm/c + Fs*Tp);                % 快时间采样点数
t = linspace(2*Rmin/c-Tp/2,2*Rmax/c+Tp/2,Nrn)';       % 快时间
R = 0.5*c*t;
Kr = B/Tp;
D = ([1:M]-1)*d;
%% load data
load Signal

data = Signal(:,1,1);
figure 
subplot(2,1,1)
plot(t,abs(data))
xlabel('t / s')
ylabel('Signal Magnitude')
subplot(2,1,2)
plot(t,angle(data))
xlabel('t / s')
ylabel('Signal Phase / rad')
% plot(t,abs(data)*cos(angle(data)))

%% 阵列综合
% 采用taylor加权
theta0 = pi*10/180;
sll = -30;                             % 旁瓣 -30dB
A =  acosh(10^abs(sll/20))/pi;
nbar = floor(2*A^2+0.5)+4;

I = taylorwin(M, nbar, sll).*exp(-j*2*pi*D'*sin(theta0)/lambda);
I = I/M;

signal = zeros(Nrn,PN);
for i = 1:M
    signal = signal + Signal(:,:,i)*I(i);
end
figure 
subplot(2,1,1)
plot(t,abs(signal(:,1)))
xlabel('t / s')
ylabel('Signal Magnitude')
subplot(2,1,2)
plot(t,angle(signal(:,1)))
xlabel('t / s')
ylabel('Signal Phase / rad')

% 方向图
angles = -90:0.2:90;
theta = pi*angles/180;
S = sum(taylorwin(M, nbar, sll).*exp(-j*2*pi*D'*sin(theta)/lambda));
S = abs(S);
figure 
plot(angles,20*log10(S/max(S)))
hold on
plot(angles,-3*ones(1,length(angles)))
axis([-90 90 -60 0])
title('阵列方向图')
xlabel('Angle / degree')
ylabel('Normalized Magnitude / dB')


%% 脉冲压缩，不进行脉冲积累
fr=([0:Nrn-1]'-ceil(Nrn/2))/Nrn*Fs;
H_f = exp(j*pi*(fr.^2/Kr))*ones(1,PN);
S_f = fftshift(fft(fftshift(signal)));
S_f = S_f.*H_f;
s_t = fftshift(ifft(fftshift(S_f)));
% mesh(abs(s_t))


%% 脉冲压缩，进行非相参积累
signal_g = sum(signal,2);
H_f_g = exp(j*pi*(fr.^2/Kr));
S_f_g = fftshift(fft(fftshift(signal_g)));
S_f_g = S_f_g.*H_f_g;
s_t_g = fftshift(ifft(fftshift(S_f_g)));
% s_t_g = sum(s_t,2);
figure
subplot(2,1,1)
plot(t,abs(s_t(:,1)))
xlabel('t / s')
ylabel('Signal Magnitude')
title('(a)无脉冲积累')
subplot(2,1,2)
plot(t,abs(s_t_g))
title('(b)非相参积累')
xlabel('t / s')
ylabel('Signal Magnitude')
