clc; clear; close all;

M      = 20;       % Number of Array Elements.
N      = 200;      % Number of Signal Samples.
n      = 1:N;      % Sample Index Vector.
lambda = 1;        % Incoming Signal Wavelength in (m).
d      = lambda/2; % Interelement Distance in (m).
SNR    = 20;       % Array Voltage Gain in dBs.
phi_s  = [20,70];       % Signal Direction angle in deg.

u_s = (d/lambda)*sin(phi_s*pi/180); % Normalized Spatial Frequency of the signal of interest.

% Desired Signal Definition:
s = zeros(M,N);

% Array snapshot (n=100) that contains the signal of interest (an impulse).
s(:,100) = (10^(SNR/20)*exp(-1i*2*pi*u_s*(0:M-1).'))/sqrt(M);

c_mf = exp(-1i*2*pi*u_s*(0:M-1).')/sqrt(M);

MC_Runs = 1;
y1 = zeros(1,N);
x1 = zeros(1,N);

for k=1:MC_Runs

    % Uncorrelated noise samples at each array element with a Gaussian distribution:
    w = (randn(M,N)+1i*randn(M,N))/sqrt(2);

    % The two signals are added to produce the overall array signal:
    x = s + w;

    % Output Calculation.
    y = c_mf'*x;
    y1 = y1 + y;

end

y_average = 1/MC_Runs*y1;
max(10*log10(abs(y_average).^2))

figure('NumberTitle', 'off','Name','Figure 11.10');
subplot(2,1,1);
% This plots the instantaneous power for every element (M waveforms).
plot(n,10*log10(abs(x).^2));
ylim([-20 25]);
grid on;
title('Instantaneous Signal Power at each Element (Ensemble of 20 Elements)');
xlabel('Sample Number');
ylabel('Output Power (dB)');

subplot(2,1,2);
plot(n,10*log10(abs(y_average).^2),'*-');
grid on;
ylim([-30 25]);
title('Instantaneous Signal Power at the Output of the Steering Vector Beamformer');
xlabel('Sample Number');
ylabel('Output Power (dB)');