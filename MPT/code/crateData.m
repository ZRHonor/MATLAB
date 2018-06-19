% create radar data cube
clear all; clc; close all; 
%% ��ʼ��
c=3e8;                                      % ���� (m/s)
lambda = 1.5;                               % ���� (m)
fc = c/lambda;                              % ��Ƶ (Hz)
B = 1e6;                                    % ��Ƶ���� (Hz)
Tp = 400e-6;                                % �������� (s)
PRT = 4000e-6;                              % �����ظ����� (s)
M = 16;                                     % ��Ԫ����
d = 0.8;                                    % ��Ԫ��� (m)
PN = 8;                                     % pulse number

% Ptarget = [200e3,    10 0];
Ptarget = [80e3,    10, -10;
           200e3,   10, -10;
           80e3,    20, -10 ];
Ptarget(:,2) = pi*Ptarget(:,2)/180;
%% ��������
Fs = 2*B;
Ts = 1/Fs;
Rmin = 50e3; Rmax = 250e3;
Rm = Rmax - Rmin;



Nrn = ceil(2*Fs*Rm/c + Fs*Tp);                % ��ʱ���������
t = linspace(2*Rmin/c-Tp/2,2*Rmax/c+Tp/2,Nrn)';       % ��ʱ��
Kr = B/Tp;
[TargetNumber,~]=size(Ptarget);
D = ([1:M]-1)*d;
%% �ز���������
% ����ά�ȣ� ��ʱ�䣬���壬��Ԫ
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