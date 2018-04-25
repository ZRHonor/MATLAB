% 信源只有一个时，检测精度与SNR成正比；
% 信源有多个时，信源之间相互干扰，检测精度与SNR不成正比

clear all
clc
close all
tic
%参数设定
M=10;
DOA =[5];
d = 0.5;
N = 400;
QAM = 16;

SNR = -3:3:15;
Mc = 100;
doa = zeros(length(SNR),3);

for i = 1:length(SNR)
    for j = 1:5
        % 建立信号模型
        X = signal(M, DOA*pi/180, N, SNR(i), QAM);
        [doa, angle, y] = DOAestimation(X, M, N, length(DOA));
        rmse(j, :) = doa - DOA;
    end
    RMSE(i,:) = sqrt(mean(rmse.^2));
end
plot(SNR,RMSE)
% plot(SNR,RMSE(:,1),SNR,RMSE(:,2),SNR,RMSE(:,3))
% legend('5','45','65')
