clear all
clc
close all
tic

DOA = [5 45 65]*pi/180;
SNR = 10;
d = 0.5;
sigma = 1/sqrt(SNR);		%the varience of the gassian noise
N = 400;
QAM = 16;
P = length(DOA);

for M = 4:2:16
    A=exp(-j*2*pi*0.5*[0:M-1]'*sin(DOA));
    %信源模型建立
    for k=1:P
        symbol = randi([0, QAM-1], 1, N);
    %     S(k,:) = qammod(symbol, QAM);
        S(k,:) = qammod(symbol, QAM);
    end
    X = awgn(A*S,SNR,'measured');

    S1=A(:,1)*A(:,1)';
    S3=A(:,3)*A(:,3)';
    n=normrnd(0,sigma,M,N);
    Cov_n=zeros(M,M);
    for m=1:N
        Cov_n=Cov_n+n(:,m)*n(:,m)';
    end
    Cov_n = Cov_n/N;
    S_I = Cov_n+(S1+S3);
    
    
    theta = 45*pi/180;
    a0 = exp(-j*2*pi*d*[0:M-1]'*sin(theta));
    a = a0/M;
    w= inv(S_I)*a0/(a0'*inv(S_I)*a0);
    % Rv=noise*noise'; 
    % Rv = mean(Rv(:));

    % SINR1 = 10*log10((W_MF'*conj(As_theta2)*As_theta2.'*W_MF)/(W_MF'*S_I*W_MF));

    SINR1(M/2-1) = 10*log10((a'*a0*a0'*a)/(a'*S_I*a))
    SINR2(M/2-1) = 10*log10((w'*a0*a0'*w)/(w'*S_I*w))
    % SINR1 = 

    
%     signal = S;
%     rsignal = A*signal;
% %     signal = S(2,:);
% %     rsignal = A(:,2)*signal;
%     X = awgn(A*S,SNR,'measured');
%     noise = X - rsignal;
% 
%     signalPower = signal*signal';
% %     signalPower = mean(temp(:));
% 
%     R = X*X'/N;
%     theta = 45*pi/180;
%     a0 = exp(-j*2*pi*d*[0:M-1]'*sin(theta))/M;
%     w = inv(R)*a0/(a0'*inv(R)*a0)/M;
%     a=a0;
% %     a=A(:,2);
% %     scatterplot(a0'*X);
% %     scatterplot(w'*X)
%     Rv=noise*noise'; 
%     Rv = mean(Rv(:));
%     SINR1(M/2-1) = 10*log10(abs(signalPower*(a0'*a)^2/(a0'*Rv*a0)));
%     SINR2(M/2-1) = 10*log10(abs(signalPower*(w'*a)^2/(w'*Rv*w)));   
% %     temp1 = abs(S*a'*a).^2;
    
end

M = 4:2:16;
plot(M,SINR1,M,SINR2);
legend('CBF算法','MVDR算法','Location','NorthWest')
xlabel('阵元个数')
ylabel('信干噪比/dB')