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

% Ptarget = [80e3,    10, 10];
Ptarget = [80e3,    10, -10;
           200e3,   10, -10;
           80e3,    20, -10 ];
Ptarget(:,2) = pi*Ptarget(:,2)/180;
%% 参数设置
Fs = 2e6;
Ts = 1/Fs;
Rmin = 50e3; Rmax = 250e3;
Rm = Rmax - Rmin;



Nrn = ceil(2*Fs*Rm/c + Fs*Tp);                % 快时间采样点数
t = linspace(2*Rmin/c-Tp/2,2*Rmax/c+Tp/2,Nrn)';       % 快时间
Kr = B/Tp;
[TargetNumber,~]=size(Ptarget);

%% 回波数据生成
% 三个维度： 快时间，脉冲，阵元
for i = 1:M
    Signal(:,:,i) = zeros(Nrn,8);
end

for i = 1:TargetNumber
    am_n=1;
    sigma=am_n*10^(Ptarget(i,3)/20);
    R_tm = Ptarget(i,1);
    % 距离向时延和相位延迟
    tau=2*R_tm/c;
    phase=pi*Kr*(t-ones(Nrn,1)*tau).^2-4*pi/lambda*ones(Nrn,1)*R_tm;
    
    Target_RangeWin=abs(t-ones(Nrn,1)*tau)<=Tp/2;
    echo_tm = sigma*exp(1j*phase).*Target_RangeWin;
    echo = repmat(echo_tm, 1, 8);
    doa = Ptarget(i,2);
    for k = 1:M
        Noise=1*am_n/sqrt(2)*(randn(Nrn,PN)+1j*randn(Nrn,PN));
        % 阵元之间的相位延迟
        ECHO(:,:,k) = echo.*exp(j*2*pi*sin(doa)*(k-1)*d/lambda) + Noise;
    end
    Signal = Signal + ECHO;
end
save Signal.mat Signal

%% 数字波束形成权矢量：泰勒加权
theta0 = pi*10/180;
sll = -30;                             % 旁瓣 -30dB
A =  acosh(10^abs(sll/20))/pi;
nbar = floor(2*A^2+0.5)+4;

I = taylorwin(M, nbar, sll);

signal = zeros(Nrn,PN);
for n =1:M
    signal = signal + Signal(:,:,n)*I(n)*exp(-j*2*pi*d*n*sin(theta0)/lambda);
end
imagesc(abs(signal));

%% 脉冲压缩 非相参积累
R = 0.5*c*t';
fr=([0:Nrn-1]'-ceil(Nrn/2))/Nrn*Fs; % 快时间变量对应的距离频率（Hz）
H_f = exp(j*pi*(fr.^2/Kr))*ones(1,PN);
S_f_tm = fftshift(fft(fftshift(signal)));
S_f_tm = S_f_tm.*H_f;
s_t_tm = fftshift(ifft(fftshift(S_f_tm)));

s_t_tm = fftshift(ifft(fftshift(S_f_tm)));
figure
imagesc(abs(s_t_tm));
colormap(gray)
title('压缩后')
figure
plot(abs(s_t_tm(:,3)),'k')
[pks locs]= findpeaks(abs(s_t_tm(:,3)'),R,'SortStr', 'descend');
locs(1:2)
% hold on
%% 脉冲压缩 相参积累
H_f = exp(j*pi*(fr.^2/Kr));
signal_s = sum(signal,2);
S_f_tm = fftshift(fft(fftshift(signal_s)));
S_f_tm = S_f_tm.*H_f;
s_t_tm = fftshift(ifft(fftshift(S_f_tm)));

s_t_tm = fftshift(ifft(fftshift(S_f_tm)));
result = abs(s_t_tm');
figure
plot(R,abs(s_t_tm),'r')
[pks locs]= findpeaks(result,R,'SortStr', 'descend');
locs(1:2)

%% 单脉冲核查测角
close all;
DATA = sum(Signal,2);
data = reshape(DATA(:,1,:),Nrn,M);

D = [1:M]*d;
angle = -90:1:75;
delta = zeros(1,length(angle));
for i = 1:length(angle)
    
    theta1 = pi*angle(i)/180;
    theta2 = pi*(angle(i)+15)/180;
    
    I1 = I.*exp(-j*2*pi*D'*sin(theta1)/lambda);
    u1 = max(abs(data*I1))
    
    I2 = I.*exp(-j*2*pi*D'*sin(theta2)/lambda);
    u2 = max(abs(data*I2))
    
    delta(i) = (u2-u1)/(u2+u1);

end
plot(angle+7.5,delta);
grid on;
close all
angle = -90:0.02:90;
theta = pi*angle/180;

S1=zeros(1,length(theta));
S2=zeros(1,length(theta));
for i = 1:M
    S1 = S1 + I(i)*exp(j*2*pi*d*i*sin(theta-pi/36)/lambda);
    S2 = S2 + I(i)*exp(j*2*pi*d*i*sin(theta+pi/36)/lambda);
end
delta = 20*log10(abs(S1)-abs(S2));
sum = 20*log10(abs(S1)+abs(S2));
plot(angle,delta)
figure
plot(angle,sum)
% theta1 = 0;
% theta2 = pi*(theta1+10)/180;
% theta1 = pi*theta1/180;
% 
% delta = 100;
% Sigma = 1;
% D = [1:M]*d;
% step = pi/36;
% 
% while(abs(delta/Sigma)>eps)
%     I1 = I.*exp(-j*2*pi*D'*n*sin(theta1)/lambda);
%     u1 = max(abs(data*I1));
%     
%     I2 = I.*exp(-j*2*pi*D'*n*sin(theta2)/lambda);
%     u2 = max(abs(data*I2));
%     
%     delta = u2-u1;
%     Sigma = u2+u1;
%     delta/Sigma
%     if(delta/Sigma>eps)
%         theta1 = theta1 + step;
%         theta2 = theta1 + pi/18;
%     else
%         theta1 = theta1 - step;
%         theta2 = theta1 + pi/18;
%     end
% end
