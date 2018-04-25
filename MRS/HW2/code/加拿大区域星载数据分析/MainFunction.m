%% ʵ�����ݴ������ô�ĳ�ۿ�����������״�����
clear all;clc;close all;
%% ���ݶ�ȡ
load DATA0CDdata1.mat;    %ĳ�ۿڳ���
data=double(data.');
[Nrn,Nan]=size(data);
c = 3e8;
%% ���״�ϵͳ��������
PRF=1256.98;           % Pulse Reputation Frequency (Hz)��ʱ�������
PRT=1/PRF;             % �����ظ����ڣ�s��
Fs=32.317e+6;          % Radar sampling rate (Hz)��ʱ��Ƶ��
fc=5.300e+9;           % Radar center frequency (Hz)�ز�Ƶ��
c=2.9979e+8;           % Speed of light (m/s)����
R0=0.0065956*c/2;      % �������ĵ�б��
Kr=0.72135e+12;        % FM rate of radar pulse (Hz/s)��ʱ���Ƶ��
Tp=41.75e-6;           % Chirp duration (s)�������ʱ��
Doppler_centroid=-6900;% ����������ƫ�ƣ�Hz��
lambda=c/fc;           % ������m
v=6.9621e+03;          %ƽ̨���ٶȣ�m/s����
fr=([0:Nrn-1]-ceil(Nrn/2))'/Nrn*Fs;     %  ����Ƶ�ʱ�����Hz��
tm=([0:Nan-1]-ceil(Nan/2))*PRT;         % ��λ��ʱ�䣨s��
fa=([0:Nan-1]-ceil(Nan/2))/Nan*PRF;     % ��λ��ʱ���Ӧ�Ķ�����Ƶ�ʣ�Hz��
tic
figure;imagesc(abs(data));
xlabel('��λ��');ylabel('������');
title('������ѹǰ�Ľ��');
%% ��������ѹ�����Ӹý���������ԵĿ������������߶�����
% step1 ����ѹ��
H_RangeComp=exp(-j*pi*fr.^2/Kr)*ones(1,Nan);
S_f_tm_RangeComp=fftshift(fft(data,[],1),1).*H_RangeComp;
s_t_tm_RangeComp=ifft(ifftshift(S_f_tm_RangeComp,1),[],1);
figure;imagesc(abs(s_t_tm_RangeComp));
xlabel('��λ��');ylabel('������');
title('after step1');
%% ���������о����߶������ͷ�λѹ�������о����������ÿ��ǣ������߶��ľ���Ϊ��(Doppler_centroid*lambda/2)*tm

% step2 ��λ����Ҷ�任
S_t_fa = fftshift(fft(s_t_tm_RangeComp,[],2),2);
figure;imagesc(abs(S_t_fa));
xlabel('��λ��');ylabel('������');
title('after step2');
% step3 �����㶯У��
% dR = fftshift(fft((Doppler_centroid*lambda/2)*tm));
% dR1 = lambda^2*R0*fa.^2/(8*v^2)\
Ka = -2*v^2/(lambda*R0);
% dR = Doppler_centroid*lambda/2*fa./(-Ka);
% H_RCMC = exp(j*4*pi*fa.*dR/c);

Nf = 2^nextpow2(Nrn);
Ns  = 2^nextpow2(Nan);
fd_r = [-Nf/2 : (Nf/2 - 1)] * Fs / Nf;
FF = ones(Ns, 1) * fd_r;                                 % FFΪN*M�ľ���
fdc = Doppler_centroid;                                                 % doppler center
fd_a = [-Ns/2 : (Ns/2 - 1)] * PRF / Ns;
FU = fd_a.' * ones(1, Nf);
Refcorr = exp(j * pi / fc^2 / Ka * (FU.*FF).^2 + j * pi * fdc^2 / fc / Ka * FF - j * pi / fc / Ka * FU.^2 .* FF); % Range-Doppler domain

S_t_fa = S_t_fa.*Refcorr;
figure;imagesc(abs(S_t_fa));
xlabel('��λ��');ylabel('������');
title('After step3');
% s_t_tm_RCMC=ifft(ifftshift(S_t_fa.*H_,2),[],2);
% S_f_tm_RCMC = S_f_tm_RangeComp.*H_RCMC;
% s_t_tm_RCMC=ifft(ifftshift(S_f_tm_RangeComp,1),[],1);
% figure;imagesc(abs(s_t_tm_RCMC));
% xlabel('��λ��');ylabel('������');
% title('�����㶯У����Ľ��');

% ��λѹ��
% S_t_fa = fftshift(fft(s_t_tm_RangeComp,[],2),2);
Ka=2*v^2/(lambda*R0);
H_AC = ones(1,Nrn)'*exp(-j*pi*fa.^2/Ka);
S_t_fa_AC = S_t_fa.*H_AC;
s_t_tm_AC = ifft(ifftshift(S_t_fa_AC,2),[],2);
figure;imagesc(abs(s_t_tm_AC));
xlabel('��λ��');ylabel('������');
title('��λѹ����Ľ��');

toc
