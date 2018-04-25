%% 实测数据处理，加拿大某港口区域的星载雷达数据
clear all;clc;close all;
%% 数据读取
load DATA0CDdata1.mat;    %某港口场景
data=double(data.');
[Nrn,Nan]=size(data);
c = 3e8;
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
%% 距离脉冲压缩，从该结果可以明显的看出发生距离走动现象
% step1 距离压缩
H_RangeComp=exp(-j*pi*fr.^2/Kr)*ones(1,Nan);
S_f_tm_RangeComp=fftshift(fft(data,[],1),1).*H_RangeComp;
s_t_tm_RangeComp=ifft(ifftshift(S_f_tm_RangeComp,1),[],1);
figure;imagesc(abs(s_t_tm_RangeComp));
xlabel('方位向');ylabel('距离向');
title('after step1');
%% 接下来进行距离走动矫正和方位压缩，其中距离弯曲不用考虑，距离走动的距离为：(Doppler_centroid*lambda/2)*tm

% step2 方位向傅里叶变换
S_t_fa = fftshift(fft(s_t_tm_RangeComp,[],2),2);
figure;imagesc(abs(S_t_fa));
xlabel('方位向');ylabel('距离向');
title('after step2');
% step3 距离徙动校正
% dR = fftshift(fft((Doppler_centroid*lambda/2)*tm));
% dR1 = lambda^2*R0*fa.^2/(8*v^2)\
Ka = -2*v^2/(lambda*R0);
% dR = Doppler_centroid*lambda/2*fa./(-Ka);
% H_RCMC = exp(j*4*pi*fa.*dR/c);

Nf = 2^nextpow2(Nrn);
Ns  = 2^nextpow2(Nan);
fd_r = [-Nf/2 : (Nf/2 - 1)] * Fs / Nf;
FF = ones(Ns, 1) * fd_r;                                 % FF为N*M的矩阵
fdc = Doppler_centroid;                                                 % doppler center
fd_a = [-Ns/2 : (Ns/2 - 1)] * PRF / Ns;
FU = fd_a.' * ones(1, Nf);
Refcorr = exp(j * pi / fc^2 / Ka * (FU.*FF).^2 + j * pi * fdc^2 / fc / Ka * FF - j * pi / fc / Ka * FU.^2 .* FF); % Range-Doppler domain

S_t_fa = S_t_fa.*Refcorr;
figure;imagesc(abs(S_t_fa));
xlabel('方位向');ylabel('距离向');
title('After step3');
% s_t_tm_RCMC=ifft(ifftshift(S_t_fa.*H_,2),[],2);
% S_f_tm_RCMC = S_f_tm_RangeComp.*H_RCMC;
% s_t_tm_RCMC=ifft(ifftshift(S_f_tm_RangeComp,1),[],1);
% figure;imagesc(abs(s_t_tm_RCMC));
% xlabel('方位向');ylabel('距离向');
% title('距离徙动校正后的结果');

% 方位压缩
% S_t_fa = fftshift(fft(s_t_tm_RangeComp,[],2),2);
Ka=2*v^2/(lambda*R0);
H_AC = ones(1,Nrn)'*exp(-j*pi*fa.^2/Ka);
S_t_fa_AC = S_t_fa.*H_AC;
s_t_tm_AC = ifft(ifftshift(S_t_fa_AC,2),[],2);
figure;imagesc(abs(s_t_tm_AC));
xlabel('方位向');ylabel('距离向');
title('方位压缩后的结果');

toc
