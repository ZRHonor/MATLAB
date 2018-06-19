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

% Ptarget = [80e3,    10, 10];
Ptarget = [80e3,    10, -10;
           200e3,   10, -10;
           80e3,    20, -10 ];
Ptarget(:,2) = pi*Ptarget(:,2)/180;
%% ��������
Fs = 2
*B;
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
    mesh(abs(echo_tm));
    for k = 1:PN
        echo = awgn(echo_tm,snr,'measured');
%         mesh(abs(echo));
%         temp = Signal(:,k,:) + reshape(echo,[Nrn,1,M]);
        Signal(:,k,:) = Signal(:,k,:) + reshape(echo,[Nrn,1,M]);
    end
end
save Signal Signal
mesh(abs(Signal(:,:,1)))
%% �źŴ��������л���
theta0 = pi*20/180;
sll = -30;                             % �԰� -30dB
A =  acosh(10^abs(sll/20))/pi;
nbar = floor(2*A^2+0.5)+4;

I = taylorwin(M, nbar, sll).*exp(-j*2*pi*D'*sin(theta0)/lambda);
% signal = zeros(Nrn,PN);
% for i = 1:M
%     signal = signal + Signal(:,:,i)*I(i);
% end
% fr=([0:Nrn-1]'-ceil(Nrn/2))/Nrn*Fs;
% H_f = exp(j*pi*(fr.^2/Kr))*ones(1,PN);
% S_f_tm = fftshift(fft(fftshift(signal)));
% S_f_tm = S_f_tm.*H_f;
% s_t_tm = fftshift(ifft(fftshift(S_f_tm)));
% 
% mesh(abs(s_t_tm));

signal = reshape(sum(Signal,2)/PN,[Nrn M]);
mesh(abs(signal))
title('abs')
figure;
mesh(angle(signal))
title('angle')

% load data.mat
% Nrn = 8000;
% signal = reshape(sum(radar_data_cube,2),[Nrn,M]);
% mesh(abs(signal))
% title('abs')
% figure;
% mesh(angle(signal))
% title('angle')



% I = exp(-j*2*pi*D'*sin(theta0)/lambda);
% I = taylorwin(M, nbar, sll).*cos(2*pi*D'*sin(theta0)/lambda);
% I = taylorwin(M, nbar, sll);
% signal = sum(signal,2);
signal = signal*I;

figure
plot(abs(signal));

fr=([0:Nrn-1]'-ceil(Nrn/2))/Nrn*Fs;
H_f = exp(j*pi*(fr.^2/Kr))*ones(1,PN);
S_f_tm = fftshift(fft(fftshift(signal)));
S_f_tm = S_f_tm.*H_f;
s_t_tm = fftshift(ifft(fftshift(S_f_tm)));
plot(abs(s_t_tm));

%% ��������

% I_bay = sin(pi*p);
% I_bay = taylorwin(16,6,-30)';
% plot(I_bay)

angle = -90:1:90;

theta = pi*angle/180;
% Bayliss ��Ȩ�Ĳ��
p = linspace(-1,1,M);
I_bay = bayliss_n(35,7,p);
S_diff = I_bay.*exp(j*2*pi*D'*(sin(theta)-sin(pi/18))/lambda);
S_diff = abs(S_diff);

% axis([-90 90 -60 0])
% Taylor ��Ȩ�ĺͲ���
I_tay = taylorwin(M,7,-30)';
S_sum = I_tay.*exp(j*2*pi*D'*(sin(theta)-sin(pi/18))/lambda);
S_sum = abs(S_sum);
figure
plot(angle,20*log10(S_sum/max(S_sum)),angle,20*log10(S_diff/max(S_diff)));
legend('�Ͳ���','���')
% title('�Ͳ���')
axis([-90 90 -60 0])
 
