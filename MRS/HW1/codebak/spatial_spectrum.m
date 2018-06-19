function [angle, ss] = spatial_spectrum(M, theta_d, SNR, N)
% 计算线阵空间谱
%   M = 线阵阵元个数
%   theta_d = 信源角度
%   SNR = 信噪比
%   N = 采样点数
    angle = -90:1:90;
    f = 1000;
    c = 1500;
    lambda=c/f;
    d=lambda/2;
%     %阵列流型A
%     P = length(theta_d);
%     for i=1:P
%         A(:,i)=exp(-j*2*pi*d*[0:M-1]'/lambda*sin(theta_d(i)));
%     end
%     %信源模型建立
%     for k=1:P
%         % S(k,:)=sqrt(10.^(snr/10))*(sind(1:N) + rand(1,N));
%         % 16QAM调制
%         QAM = 16;
%         x = randi(QAM,1,N) - 1;
%         y = qammod(x, QAM);
%         S(k,:) = awgn(y, SNR);
%     end
%
%     %接收信号模型建立
%     X=A*S;
    X = signal(M, theta_d, N, SNR);
%     a = exp(-j*2*pi*d*[0:M-1]'/lambda*sin(5));
%     sig = a'*X;
%     scatterplot(sig)
    %协方差矩阵
    R=X*X'/N;

    syms theta
    a = exp(-j*2*pi*d*[0:M-1]'/lambda*sin(theta));
    y(theta) = sqrt(abs(a'*R*a));
    y_theta = real(abs(subs(y, theta, angle*pi/180)));
    ss = real(double(y_theta));
end