%------���ܣ�SAR��Ŀ����棬Ŀ���б�����̽�����б��-��λ��άƽ��Ļ�����-----------
%------Ŀ�����Ĳ����������������ѹ���������㶯�����Լ���λѹ��
clear all; clc; close all;
%% �������״�ϵͳ��������
c=3e8;                                     % ���٣�m/s��
fc=8.85e9;                                 % ��Ƶ��Hz��
lambda=c/fc;                               % ������m��
H=6000;                                    % ƽ̨�߶ȣ�m��
R0=10000;                                  % Ŀ�����б��
v=120;                                     % ƽ̨�ٶȣ�m/s��
%% ������Ŀ�����
thetaAz=0;                                 % ��λ��Ƕ�
D=1;
theta_azi=lambda/D;                        % ��λ�������
AL=R0*theta_azi ;                          % ���ߺϳɿ׾�����
Ta=AL/v;                                   % ��λ�ϳɿ׾�ʱ��
PRF=1000;                                  % �����ظ�Ƶ��(Hz)
PRT=1/PRF;                                 % �����ظ�����(s)
Nan=ceil(1.5*AL/v/PRT);                    % ��ʱ���������
%----��ʱ��˵������������AL/v/PRT������Ŀ������Ȳ��ںϳɿ׾�ʱ���ڣ�
%----�������з�λ�Ӵ�������˻����ܹ���֤Ŀ����һ����Ч�ĺϳɿ׾�ʱ����
tm=([0:Nan-1]-ceil(Nan/2))*PRT;            % ��ʱ�����(s)
fa=([0:Nan-1]-ceil(Nan/2))/Nan*PRF;        % ��ʱ���Ӧ�Ķ�����Ƶ�ʱ���(Hz)
%----------��ʱ����----------%
Tp=0.2e-6;                                 % ������(s)
B=100e6;                                   % �źŴ���(Hz)
Kr=B/Tp;                                   % �źŵ�Ƶ��(Hz/s)
Fs=B*1.2;                                  % ��ʱ������ʣ�Hz��
dt=1/Fs;                                   % ��ʱ����������s��
Rm=100;                                    % ����ά��Զб���ȥ���б�ࣨm��
Nrn=ceil(Rm/(c/2/Fs)+Tp/dt);               % ��ʱ���������
t=linspace(2*(R0-Rm/2)/c-Tp/2,2*(R0+Rm/2)/c+Tp/2,Nrn).';% ��ʱ���������ʱ�䣨s��
%- R0������άб������б�࣬Rm����Զб���ȥ���б��
fr=([0:Nrn-1]'-ceil(Nrn/2))/Nrn*Fs;        % ��ʱ�������Ӧ�ľ���Ƶ�ʣ�Hz��
%% Ŀ���˶��������ã����ɸĶ���
Range_Azi_initialize=0;
%             [��λ�ؾࣨm����        ���б�ࣨm����   SNR��dB��]
Ptarget =[Range_Azi_initialize,     R0+c/2/B*0,         15;
          Range_Azi_initialize,     R0+c/2/B*10,        15;
         Range_Azi_initialize+50,   R0+c/2/B*15,        15;
         Range_Azi_initialize-50,   R0+c/2/B*20,        15];
%% �źŲ���
[TargetNumber,~]=size(Ptarget);
Echo=zeros(Nrn,Nan);
for ii=1:TargetNumber
    am_n=1;
    sigma=am_n*10^(Ptarget(ii,3)/20);
    x0_temp=Ptarget(ii,1);
    R0_temp=Ptarget(ii,2);
    disp(['Ŀ��',num2str(ii),'�ľ���������Ԫ��Ϊ��',num2str(R0_temp*lambda^2/32/(D/2)^2/(c/2/Fs))])
    R_tm=sqrt(R0_temp^2+(v*tm-x0_temp).^2);
    tau=2*R_tm/c;
    phase=pi*Kr*(t*ones(1,Nan)-ones(Nrn,1)*tau).^2-4*pi/lambda*ones(Nrn,1)*R_tm;
    Target_RangeWin=abs(t*ones(1,Nan)-ones(Nrn,1)*tau)<=Tp/2;
    Target_AziWin=abs((v*ones(Nrn,1)*tm)-Ptarget(ii,1))<=AL/2;
    Echo=Echo+sigma*exp(1j*phase).*Target_RangeWin.*Target_AziWin;
end
Noise=1*am_n/sqrt(2)*(randn(Nrn,Nan)+1j*randn(Nrn,Nan));  %---����������
Signal=Echo+Noise;                                       %---������Ŀ��+�����ز��ź�
figure;imagesc(abs(Signal));
%colorbar;
xlabel('��λ��');ylabel('������');
title('�ز�����');
colormap(gray)

%% ���ݴ���
tic

% step1 ����ѹ��
S_f_tm = fftshift(fft(fftshift(Signal)));
H_f = exp(j*pi*(fr.^2/Kr))*ones(1,Nan);
S_f_tm = S_f_tm.*H_f;
s_t_tm = fftshift(ifft(fftshift(S_f_tm)));

figure;imagesc(abs(s_t_tm));
%colorbar;
xlabel('��λ��');ylabel('������');
title('�����㶯У��ǰ�Ľ��');
% title('a')
colormap(gray);


% step2 ��λFFT
S_t_fa = fftshift(fft(fftshift(s_t_tm, 2), [], 2), 2);

% step3 �����㶯У��
S_f_fa= fftshift(fft(fftshift(S_t_fa)));

R = linspace(R0-Rm/2-Tp*c/4,R0+Rm/2+Tp*c/4,Nrn)';
dR_fa = (lambda^2/(8*v^2))*(R)*(fa.^2);
S_f_fa = S_f_fa.*exp(j*4*pi*(fr*ones(1,Nan)).*dR_fa/c);

S_t_fa = fftshift(ifft(fftshift(S_f_fa)));

s_t_tm = fftshift(ifft(S_t_fa, [], 2), 2);
figure;imagesc(abs(s_t_tm));
%colorbar;
xlabel('��λ��');ylabel('������');
title('�����㶯У��֮��Ľ��');
% title('b')
colormap(gray);


% step4 ��λѹ��
Ka=2*v^2/(lambda*R0);
Haz_t_fa = ones(1,Nrn)'*exp(-j*pi*fa.^2/Ka);
S_t_fa = S_t_fa.*Haz_t_fa;

% step5 ��λIFFT
s_t_tm = fftshift(ifft(S_t_fa, [], 2), 2);
figure;imagesc(abs(s_t_tm));
%colorbar;
xlabel('��λ��');ylabel('������');
title('Result');
colormap(gray);

toc