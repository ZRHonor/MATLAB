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
Noise=0*am_n/sqrt(2)*(randn(Nrn,Nan)+1j*randn(Nrn,Nan));  %---����������
Signal=Echo;                                       %---������Ŀ��+�����ز��ź�
figure;imagesc(abs(Signal));colorbar;
xlabel('��λ��');ylabel('������');
title('������ѹǰ������');
%% ���������о�������ѹ���������㶯�����ͷ�λѹ���Ȳ���

%% edit on 4.14
tic

x = R0 + linspace(-Rm/2, Rm/2, Nrn);
phi0 = x'*(1./sqrt(fc^2 - fa.^2));
phi1 = fc*x'*(1./sqrt(fc^2 - fa.^2));
phi2 = 1/2*x'*(fa.^2./(fc^2 - fa.^2).^1.5);
Cs = ones(1, Nrn)'*(fc./sqrt(fc^2 - fa.^2 ) - 1);
Ks = 1./(1/Kr  - 2*phi2);

% Step1: ��λ��FFT
S_t_fa = fftshift(fft(fftshift(Signal, 2), [], 2), 2);
figure;imagesc(abs(S_t_fa));colorbar;
xlabel('fa');ylabel('t');
title('step1');
% Step2: Chirp Scaling
Scs_t_fa = S_t_fa.*exp(1j*pi*Cs.*Ks.*(x'*ones(1,Nan) - R0*(1 + Cs).^2));
figure;imagesc(abs(S_t_fa));colorbar;
xlabel('fa');ylabel('t');
title('step2');

% Step3: ������FFT
Scs_fr_fa = fftshift(fft(fftshift(Scs_t_fa)));
Scs_t_fa = S_t_fa.*exp(j*pi*Cs.*Ks.*(x'*ones(1,Nan) - R0*(1 + Cs).^2));
figure;imagesc(abs(Scs_t_fa));colorbar;
xlabel('fa');ylabel('fr');
title('step3');

% Step4: ������Ǩ��У���;�����ƥ���˲�
Srmc_fr_fa = Scs_fr_fa.*exp(j*pi*(fr'.^2'*ones(1,Nan))./(1+Cs)./Ks+j*2*pi*R0*Cs.*(fr''*ones(1,Nan)));
figure;imagesc(abs(Srmc_fr_fa));colorbar;
xlabel('fa');ylabel('fr');
title('step4');

% Step5: ������IFFT
Srmc_t_fa = fftshift(ifft(fftshift(Srmc_fr_fa)));
figure;imagesc(abs(Srmc_t_fa));colorbar;
xlabel('fa');ylabel('t');
title('step5');

% Step6: ��λ��ƥ���˲�
Ka=2*v^2/(lambda*R0);
Haz_t_fa = ones(1,Nrn)'*exp(-j*pi*fa.^2/Ka);
Srmc_t_fa = Srmc_t_fa.*Haz_t_fa;
% Srmc_t_fa = Srmc_t_fa.*exp(-j*pi*Ks.*Cs.*(1+Cs).*((x-R0).^2'*ones(1,Nan))-j*2*pi*phi0);
figure;imagesc(abs(Srmc_t_fa));colorbar;
xlabel('fa');ylabel('t');
title('step6');

% Step7: ��λ��IFFT
% s_t_ta = ifftshift(ifft(ifftshift(Srmc_t_fa, 2), [], 2)) ;
s_t_ta = fftshift(ifft(Srmc_t_fa, [], 2), 2);
figure;imagesc(abs(s_t_ta));colorbar;
xlabel('��λ��');ylabel('������');
title('result');
toc