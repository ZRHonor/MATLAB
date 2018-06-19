%% 实测数据处理，加拿大某港口区域的星载雷达数据
clear all;clc;close all;
%% 数据读取
load DATA0CDdata1.mat;    %某港口场景
data=double(data.');
[Nrn,Nan]=size(data);
% c = 3e8;
%% 该雷达系统参数配置
PRF=1256.98;           % Pulse Reputation Frequency (Hz)慢时间采样率
PRT=1/PRF;             % 脉冲重复周期（s）
Fs=32.317e+6;          % Radar sampling rate (Hz)快时间频率
fc=5.300e+9;           % Radar center frequency (Hz)载波频率
c=2.9979e+8;           % Speed of light (m/s)光速
R0=0.0065956*c/2;      % 场景中心的斜距
Kr=0.72135e+12;        % FM rate of radar pulse (Hz/s)快时间调频率
Tp=41.75e-6;           % Chirp duration (s)脉冲持续时间
Doppler_centroid=-6900;% 多普勒中心偏移（Hz）
lambda=c/fc;           % 波长（m
v=6.9621e+03;          %平台的速度（m/s））
fr=([0:Nrn-1]-ceil(Nrn/2))'/Nrn*Fs;     %  距离频率变量（Hz）
tm=([0:Nan-1]-ceil(Nan/2))*PRT;         % 方位慢时间（s）
fa=([0:Nan-1]-ceil(Nan/2))/Nan*PRF;     % 方位慢时间对应的多普勒频率（Hz）
tic
figure;imagesc(abs(data));
xlabel('方位向');ylabel('距离向');
title('距离脉压前的结果');


% step1 距离压缩
H_RangeComp=exp(-j*pi*fr.^2/Kr)*ones(1,Nan);
S_f_tm=fftshift(fft(data,[],1),1).*H_RangeComp;

beforeRCMC=ifft(ifftshift(S_f_tm,1),[],1);
figure;imagesc(abs(beforeRCMC));
title('距离徙动校正前')
colormap(gray)
% title('a')
xlabel('方位向');ylabel('距离向');

% step2 距离徙动校正
dR = (Doppler_centroid*lambda/2)*tm;
S_f_tm = S_f_tm.*exp(-j*4*pi*fr*dR/c);

s_t_tm=ifft(ifftshift(S_f_tm,1),[],1);
figure;imagesc(abs(s_t_tm));
colormap(gray)
title('距离徙动校正后')
% title('b')
xlabel('方位向');ylabel('距离向');

% step3 方位向傅里叶变换
S_t_fa = fft(fftshift(s_t_tm,2),[],2);

% step4 方位压缩
Rm = Nrn/Fs*c;
R = linspace(R0-Rm/2,R0+Rm/2,Nrn)';
Ka=2*v^2./(lambda*R);
H_AC = exp(-j*pi*(1./Ka)*fa.^2);
% Ka=2*v^2./(lambda*R0);
% H_AC = ones(1,Nrn)'*exp(-j*pi*(1./Ka)*fa.^2);
S_t_fa_AC = S_t_fa.*H_AC;
s_t_tm_AC = fftshift(ifft(fftshift(S_t_fa_AC,2),[],2),2);

result = abs(s_t_tm_AC);
% result(result>50)=50;
figure;imagesc(result);
axis equal;axis off;
colormap(gray);

imcontrast;
xlabel('方位向');ylabel('距离向');
title('result');

toc
