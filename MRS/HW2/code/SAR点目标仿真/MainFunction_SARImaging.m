%------功能：SAR点目标仿真，目标的斜距历程建立在斜距-方位二维平面的基础上-----------
%------目标成像的步骤包括：距离脉冲压缩，距离徙动矫正以及方位压缩
clear all; clc; close all;
%% 基本的雷达系统参数设置
c=3e8;                                     % 光速（m/s）
fc=8.85e9;                                 % 载频（Hz）
lambda=c/fc;                               % 波长（m）
H=6000;                                    % 平台高度（m）
R0=10000;                                  % 目标最近斜距
v=120;                                     % 平台速度（m/s）
%% 基本的目标参数
thetaAz=0;                                 % 方位向角度
D=1;
theta_azi=lambda/D;                        % 方位向波束宽度
AL=R0*theta_azi ;                          % 天线合成孔径长度
Ta=AL/v;                                   % 方位合成孔径时间
PRF=1000;                                  % 脉冲重复频率(Hz)
PRT=1/PRF;                                 % 脉冲重复周期(s)
Nan=ceil(1.5*AL/v/PRT);                    % 慢时间采样个数
%----慢时间说明：理论上是AL/v/PRT，但是目标可能先不在合成孔径时间内，
%----但后面有方位加窗处理，因此还是能够保证目标在一个有效的合成孔径时间内
tm=([0:Nan-1]-ceil(Nan/2))*PRT;            % 慢时间变量(s)
fa=([0:Nan-1]-ceil(Nan/2))/Nan*PRF;        % 慢时间对应的多普勒频率变量(Hz)
%----------快时间域----------%
Tp=0.2e-6;                                 % 脉冲宽度(s)
B=100e6;                                   % 信号带宽(Hz)
Kr=B/Tp;                                   % 信号调频率(Hz/s)
Fs=B*1.2;                                  % 快时间采样率（Hz）
dt=1/Fs;                                   % 快时间采样间隔（s）
Rm=100;                                    % 俯仰维最远斜距减去最近斜距（m）
Nrn=ceil(Rm/(c/2/Fs)+Tp/dt);               % 快时间采样个数
t=linspace(2*(R0-Rm/2)/c-Tp/2,2*(R0+Rm/2)/c+Tp/2,Nrn).';% 快时间采样数据时间（s）
%- R0：俯仰维斜距中心斜距，Rm：最远斜距减去最近斜距
fr=([0:Nrn-1]'-ceil(Nrn/2))/Nrn*Fs;        % 快时间变量对应的距离频率（Hz）
%% 目标运动参数设置（均可改动）
Range_Azi_initialize=0;
%             [方位地距（m），        最近斜距（m），   SNR（dB）]
Ptarget =[Range_Azi_initialize,     R0+c/2/B*0,         15;
          Range_Azi_initialize,     R0+c/2/B*10,        15;
         Range_Azi_initialize+50,   R0+c/2/B*15,        15;
         Range_Azi_initialize-50,   R0+c/2/B*20,        15];
%% 信号产生
[TargetNumber,~]=size(Ptarget);
Echo=zeros(Nrn,Nan);
for ii=1:TargetNumber
    am_n=1;
    sigma=am_n*10^(Ptarget(ii,3)/20);
    x0_temp=Ptarget(ii,1);
    R0_temp=Ptarget(ii,2);
    disp(['目标',num2str(ii),'的距离弯曲单元数为：',num2str(R0_temp*lambda^2/32/(D/2)^2/(c/2/Fs))])
    R_tm=sqrt(R0_temp^2+(v*tm-x0_temp).^2);
    tau=2*R_tm/c;
    phase=pi*Kr*(t*ones(1,Nan)-ones(Nrn,1)*tau).^2-4*pi/lambda*ones(Nrn,1)*R_tm;
    Target_RangeWin=abs(t*ones(1,Nan)-ones(Nrn,1)*tau)<=Tp/2;
    Target_AziWin=abs((v*ones(Nrn,1)*tm)-Ptarget(ii,1))<=AL/2;
    Echo=Echo+sigma*exp(1j*phase).*Target_RangeWin.*Target_AziWin;
end
Noise=0*am_n/sqrt(2)*(randn(Nrn,Nan)+1j*randn(Nrn,Nan));  %---产生的噪声
Signal=Echo;                                       %---产生的目标+噪声回波信号
figure;imagesc(abs(Signal));colorbar;
xlabel('方位向');ylabel('距离向');
title('距离脉压前的数据');
%% 接下来进行距离脉冲压缩，距离徙动矫正和方位压缩等步骤

%% edit on 4.14
tic

x = R0 + linspace(-Rm/2, Rm/2, Nrn);
phi0 = x'*(1./sqrt(fc^2 - fa.^2));
phi1 = fc*x'*(1./sqrt(fc^2 - fa.^2));
phi2 = 1/2*x'*(fa.^2./(fc^2 - fa.^2).^1.5);
Cs = ones(1, Nrn)'*(fc./sqrt(fc^2 - fa.^2 ) - 1);
Ks = 1./(1/Kr  - 2*phi2);

% Step1: 方位向FFT
S_t_fa = fftshift(fft(fftshift(Signal, 2), [], 2), 2);
figure;imagesc(abs(S_t_fa));colorbar;
xlabel('fa');ylabel('t');
title('step1');
% Step2: Chirp Scaling
Scs_t_fa = S_t_fa.*exp(1j*pi*Cs.*Ks.*(x'*ones(1,Nan) - R0*(1 + Cs).^2));
figure;imagesc(abs(S_t_fa));colorbar;
xlabel('fa');ylabel('t');
title('step2');

% Step3: 距离向FFT
Scs_fr_fa = fftshift(fft(fftshift(Scs_t_fa)));
Scs_t_fa = S_t_fa.*exp(j*pi*Cs.*Ks.*(x'*ones(1,Nan) - R0*(1 + Cs).^2));
figure;imagesc(abs(Scs_t_fa));colorbar;
xlabel('fa');ylabel('fr');
title('step3');

% Step4: 距离向迁移校正和距离向匹配滤波
Srmc_fr_fa = Scs_fr_fa.*exp(j*pi*(fr'.^2'*ones(1,Nan))./(1+Cs)./Ks+j*2*pi*R0*Cs.*(fr''*ones(1,Nan)));
figure;imagesc(abs(Srmc_fr_fa));colorbar;
xlabel('fa');ylabel('fr');
title('step4');

% Step5: 距离向IFFT
Srmc_t_fa = fftshift(ifft(fftshift(Srmc_fr_fa)));
figure;imagesc(abs(Srmc_t_fa));colorbar;
xlabel('fa');ylabel('t');
title('step5');

% Step6: 方位向匹配滤波
Ka=2*v^2/(lambda*R0);
Haz_t_fa = ones(1,Nrn)'*exp(-j*pi*fa.^2/Ka);
Srmc_t_fa = Srmc_t_fa.*Haz_t_fa;
% Srmc_t_fa = Srmc_t_fa.*exp(-j*pi*Ks.*Cs.*(1+Cs).*((x-R0).^2'*ones(1,Nan))-j*2*pi*phi0);
figure;imagesc(abs(Srmc_t_fa));colorbar;
xlabel('fa');ylabel('t');
title('step6');

% Step7: 方位向IFFT
% s_t_ta = ifftshift(ifft(ifftshift(Srmc_t_fa, 2), [], 2)) ;
s_t_ta = fftshift(ifft(Srmc_t_fa, [], 2), 2);
figure;imagesc(abs(s_t_ta));colorbar;
xlabel('方位向');ylabel('距离向');
title('result');
toc