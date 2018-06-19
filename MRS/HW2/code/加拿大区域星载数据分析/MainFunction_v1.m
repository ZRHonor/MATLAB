%% ʵ�����ݴ������ô�ĳ�ۿ�����������״�����
clear all;clc;close all;
%% ���ݶ�ȡ
load DATA0CDdata1.mat;    %ĳ�ۿڳ���
data=double(data.');
[Nrn,Nan]=size(data);
% c = 3e8;
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


% step1 ����ѹ��
H_RangeComp=exp(-j*pi*fr.^2/Kr)*ones(1,Nan);
S_f_tm=fftshift(fft(data,[],1),1).*H_RangeComp;

beforeRCMC=ifft(ifftshift(S_f_tm,1),[],1);
figure;imagesc(abs(beforeRCMC));
title('�����㶯У��ǰ')
colormap(gray)
% title('a')
xlabel('��λ��');ylabel('������');

% step2 �����㶯У��
dR = (Doppler_centroid*lambda/2)*tm;
S_f_tm = S_f_tm.*exp(-j*4*pi*fr*dR/c);

s_t_tm=ifft(ifftshift(S_f_tm,1),[],1);
figure;imagesc(abs(s_t_tm));
colormap(gray)
title('�����㶯У����')
% title('b')
xlabel('��λ��');ylabel('������');

% step3 ��λ����Ҷ�任
S_t_fa = fft(fftshift(s_t_tm,2),[],2);

% step4 ��λѹ��
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
xlabel('��λ��');ylabel('������');
title('result');

toc
