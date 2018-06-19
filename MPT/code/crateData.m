% create radar data cube
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

% Ptarget = [200e3,    10 0];
Ptarget = [80e3,    10, -10;
           200e3,   10, -10;
           80e3,    20, -10 ];
Ptarget(:,2) = pi*Ptarget(:,2)/180;
%% 参数设置
Fs = 2*B;
Ts = 1/Fs;
Rmin = 50e3; Rmax = 250e3;
Rm = Rmax - Rmin;



Nrn = ceil(2*Fs*Rm/c + Fs*Tp);                % 快时间采样点数
t = linspace(2*Rmin/c-Tp/2,2*Rmax/c+Tp/2,Nrn)';       % 快时间
Kr = B/Tp;
[TargetNumber,~]=size(Ptarget);
D = ([1:M]-1)*d;
%% 回波数据生成
% 三个维度： 快时间，脉冲，阵元
for i = 1:M
    Signal(:,:,i) = zeros(Nrn,8);
end


for i = 1:TargetNumber
    doa = Ptarget(i,2);
    snr = Ptarget(i,3);
    R_tm = 2*Ptarget(i,1) + D*sin(doa);
    tau = R_tm/c;
    temp = R_tm-floor(R_tm/lambda)*lambda;
    phase = pi*Kr*(t-ones(Nrn,1)*tau).^2 + 2*pi*R_tm/lambda;
    Target_RangeWin=abs(t-ones(Nrn,1)*tau)<=Tp/2;
    echo_tm = exp(j*phase).*Target_RangeWin;
    for k = 1:PN
        echo = awgn(echo_tm,snr,'measured');
        Signal(:,k,:) = Signal(:,k,:) + reshape(echo,[Nrn,1,M]);
    end
end
save Signal Signal