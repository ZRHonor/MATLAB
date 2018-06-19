clear all; clc;
c=3e8;                                      % ���� (m/s)
lambda = 1.5;                               % ���� (m)
fc = c/lambda;                              % ��Ƶ (Hz)
B = 1e6;                                    % ��Ƶ���� (Hz)
Tp = 400e-6;                                % �������� (s)
PRT = 4000e-6;                              % �����ظ����� (s)
M = 16;                                     % ��Ԫ����
d = 0.8;                                    % ��Ԫ��� (m)
PN = 8;                                     % pulse number

%% ��������
Fs = 2*B;
Ts = 1/Fs;
Rmin = 50e3; Rmax = 250e3;
Rm = Rmax - Rmin;

Nrn = ceil(2*Fs*Rm/c + Fs*Tp);                % ��ʱ���������
t = linspace(2*Rmin/c-Tp/2,2*Rmax/c+Tp/2,Nrn)';       % ��ʱ��
R = 0.5*c*t;
Kr = B/Tp;
D = ([1:M]-1)*d;


%% ��������
angles = -90:0.2:90;
theta = pi*angles/180;
% �Ͳ���������taylor��Ȩ
sll = -30;
A =  acosh(10^abs(sll/20))/pi;
nbar = floor(2*A^2+0.5)+4;
I_tay = taylorwin(M, nbar, sll);
S_sum = sum(I_tay.*exp(-j*2*pi*D'*(sin(theta))/lambda));

% ���������bayliss��Ȩ
p = linspace(-1,1,M);
I_bay = bayliss_n(30,7,p)';
S_diff = sum(I_bay.*exp(-j*2*pi*D'*(sin(theta))/lambda));

% ��ͼ
S_diff = abs(S_diff);
S_sum = abs(S_sum);
figure
plot(angles,20*log10(S_diff/max(S_diff)),angles,20*log10(S_sum/max(S_sum)));
axis([-90 90 -60 0])
legend('���','�Ͳ���')
xlabel('Angle / degree')
ylabel('Normalized Magnitude / dB')


% ��������
angles = -1.5:0.01:1.5;
thetas = pi*angles/180;
Sigma = zeros(1,length(angles));
Delta = zeros(1,length(angles));
for i = 1:length(angles)
    data = exp(j*2*pi*D*sin(thetas(i))/lambda);
    I = I_tay;
    Sigma(i) = data*I;
    I = I_bay;
    Delta(i) = data*I;
end
% ��������
ratio = imag(Delta./Sigma);
% ��ϼ�������
b=polyfit(angles,ratio,1);
km=b(1);
yy = polyval(b,angles);
figure
plot(angles,ratio,angles,yy)
legend('��������','��ϼ�������','Location','NorthEast')
xlabel('\theta_{\it{t}} / degree')
ylabel('\Delta/\Sigma')

%% ����ʱ�̵��źŲ��
load Signal
signal = reshape(Signal(:,1,:),[Nrn M]);
data = signal(2401,:);                      % ȡһ��ʱ�̵��ź�

% �Ͳ���ɨ�裬Ѱ����Ӧ���ֵ����ΪԤ��Ŀ��Ƕȣ������ź���
angles = 5:0.001:15;
thetas = pi*angles/180;
Sigma = zeros(1,length(angles));
Delta = zeros(1,length(angles));
for i = 1:length(angles)
    I = I_tay.*exp(-j*2*pi*D'*(sin(thetas(i)))/lambda);
    Sigma(i) = data*I;
    I = I_bay.*exp(-j*2*pi*D'*(sin(thetas(i)))/lambda);
    Delta(i) = data*I;
end
figure 
plot(angles,abs(Sigma))
xlabel('\theta / degree')
ylabel('Response Magnitude')
% plotyy(angles,abs(Sigma),angles,imag(Sigma))
figure
plot(angles,abs(Delta))
xlabel('\theta / degree')
ylabel('Response Magnitude')
temp = (abs(Sigma) == max(abs(Sigma)));
Sigma0 = Sigma(temp);
theta0 = thetas(temp);
angle0 = angles(temp);

I = I_bay.*exp(-j*2*pi*D'*(sin(theta0))/lambda);
Delta0 = data*I;

angle_t = imag(Delta0/Sigma0)/km;
disp('����ʱ�̵Ĳ�ǽ��')
angle_fin = angle0 - angle_t

%% ����5��ʱ�̵��źŽ��в��,
result = zeros(1,9);
for k = 1:9
    data = signal(2399+k,:);                      % ȡһ��ʱ�̵��ź�

    % �Ͳ���ɨ�裬Ѱ����Ӧ���ֵ����ΪԤ��Ŀ��Ƕȣ������ź���
    angles = 5:0.001:15;
    thetas = pi*angles/180;
    Sigma = zeros(1,length(angles));
    Delta = zeros(1,length(angles));
    for i = 1:length(angles)
        I = I_tay.*exp(-j*2*pi*D'*(sin(thetas(i)))/lambda);
        Sigma(i) = data*I;
        I = I_bay.*exp(-j*2*pi*D'*(sin(thetas(i)))/lambda);
        Delta(i) = data*I;
    end
    temp = (abs(Sigma) == max(abs(Sigma)));
    Sigma0 = Sigma(temp);
    theta0 = thetas(temp);
    angle0 = angles(temp);

    I = I_bay.*exp(-j*2*pi*D'*(sin(theta0))/lambda);
    Delta0 = data*I;

    angle_t = imag(Delta0/Sigma0)/km;
    result(k) = angle0 - angle_t;
end
disp('���ʱ��ƽ���Ĳ�ǽ��')
mean(result)
