% 信源只有一个时，检测精度与SNR成正比；
% 信源有多个时，信源之间相互干扰，检测精度与SNR不成正比

clear all
clc
close all
tic
%参数设定
M=10;
DOA =[45];
d = 0.5;
N = 400;
QAM = 16;

SNR = -3:3:15;
Mc = 100;
% doa = zeros(length(SNR),3);
timer=0;

for i = 1:length(SNR)
    for j = 1:100
        X = signal(M, DOA*pi/180, N, SNR(i), QAM);
        R = X*X'/N;
        [doa_CBF, angle] = DOAestimation(X, M, N, length(DOA),R);
        doa_CBF
        rmse_CBF(j, :) = doa_CBF - DOA;
%         rmse_MVDR(j, :) = doa_MVDR - DOA;
        timer=timer+1
    end
    RMSE_CBF(i,:) = sqrt(mean(rmse_CBF.^2));
%     RMSE_MVDR(i,:) = sqrt(mean(rmse_MVDR.^2));
end
% plot(SNR,RMSE_CBF,SNR,RMSE_MVDR)
save RMSE_CBF
% save RMSE_MVDR
% plot(SNR,RMSE_CBF(:,1),SNR,RMSE_CBF(:,2),SNR,RMSE_CBF(:,3), ...
%      SNR,RMSE_MVDR(:,1),SNR,RMSE_MVDR(:,2),SNR,RMSE_MVDR(:,3))
% legend('CBF-5°','CBF-45°','CBF-65°','MVDR-5°','MVDR-45°','MVDR-65°')
plot(SNR,RMSE_CBF);
xlabel('SNR/dB')
ylabel('RMSE')
