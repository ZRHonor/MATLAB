% ��Դֻ��һ��ʱ����⾫����SNR�����ȣ�
% ��Դ�ж��ʱ����Դ֮���໥���ţ���⾫����SNR��������

clear all
clc
close all
tic
%�����趨
M=10;
DOA =[45];
d = 0.5;
N = 400;
QAM = 16;

SNR = -3:3:15;
Mc = 100;

for i = 1:length(SNR)
    for j = 1:100
        X = signal(M, DOA*pi/180, N, SNR(i), QAM);
        [doa_CBF, angle] = DOAestimation(X, M, N, length(DOA));
        rmse_CBF(j, :) = doa_CBF - DOA;
    end
    RMSE_CBF(i,:) = sqrt(mean(rmse_CBF.^2));
end

save RMSE_CBF
plot(SNR,RMSE_CBF);
xlabel('SNR/dB')
ylabel('RMSE')

toc
