clear all
clc
close all
tic
%参数设定
M = 10;
DOA = [5 45 65]*pi/180;
SNR = 10;
sigma = 1/sqrt(SNR);		%the varience of the gassian noise
d = 0.5;
N = 400;
QAM = 16;
P = length(DOA);
A=exp(-j*2*pi*0.5*[0:M-1]'*sin(DOA));
%信源模型建立
for k=1:P
    symbol = randi([0, QAM-1], 1, N);
%     S(k,:) = symbol;
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



% signal = S(2,:);
% 
% noise = X - A(:,2)*signal;
% 
% signalPower = signal*signal';
% % signalPower = mean(temp(:));
% 
% 
% 
% 
% R = X*X'/N;
theta = 45*pi/180;
a0 = exp(-j*2*pi*d*[0:M-1]'*sin(theta));
a = a0/M;
w= inv(S_I)*a0/(a0'*inv(S_I)*a0);
% Rv=noise*noise'; 
% Rv = mean(Rv(:));

% SINR1 = 10*log10((W_MF'*conj(As_theta2)*As_theta2.'*W_MF)/(W_MF'*S_I*W_MF));

SINR1 = 10*log10((a'*a0*a0'*a)/(a'*S_I*a))
SINR2 = 10*log10((w'*a0*a0'*w)/(w'*S_I*w))
% SINR1 = 
% SINR1 = real(10*log10(signalPower*abs(a'*a)^2/(a'*Rv*a)))
% SINR2 = real(10*log10(signalPower*abs(w'*a)^2/(w'*Rv*w)))
