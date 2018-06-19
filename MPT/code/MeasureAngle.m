clear all; clc;
c=3e8;                                      % 光速 (m/s)
lambda = 1.5;                               % 波长 (m)
fc = c/lambda;                              % 载频 (Hz)
B = 1e6;                                    % 调频带宽 (Hz)
Tp = 400e-6;                                % 发射脉宽 (s)
PRT = 4000e-6;                              % 脉冲重复周期 (s)
M = 16;                                     % 阵元个数
d = 0.8;                                    % 阵元间距 (m)
PN = 8;                                     % pulse number

%% 参数设置
Fs = 2*B;
Ts = 1/Fs;
Rmin = 50e3; Rmax = 250e3;
Rm = Rmax - Rmin;

Nrn = ceil(2*Fs*Rm/c + Fs*Tp);                % 快时间采样点数
t = linspace(2*Rmin/c-Tp/2,2*Rmax/c+Tp/2,Nrn)';       % 快时间
R = 0.5*c*t;
Kr = B/Tp;
D = ([1:M]-1)*d;


%% 单脉冲测角
angles = -90:0.2:90;
theta = pi*angles/180;
% 和波束，采用taylor加权
sll = -30;
A =  acosh(10^abs(sll/20))/pi;
nbar = floor(2*A^2+0.5)+4;
I_tay = taylorwin(M, nbar, sll);
S_sum = sum(I_tay.*exp(-j*2*pi*D'*(sin(theta))/lambda));

% 差波束，采用bayliss加权
p = linspace(-1,1,M);
I_bay = bayliss_n(30,7,p)';
S_diff = sum(I_bay.*exp(-j*2*pi*D'*(sin(theta))/lambda));

% 画图
S_diff = abs(S_diff);
S_sum = abs(S_sum);
figure
plot(angles,20*log10(S_diff/max(S_diff)),angles,20*log10(S_sum/max(S_sum)));
axis([-90 90 -60 0])
legend('差波束','和波束')
xlabel('Angle / degree')
ylabel('Normalized Magnitude / dB')


% 鉴角曲线
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
% 鉴角曲线
ratio = imag(Delta./Sigma);
% 拟合鉴角曲线
b=polyfit(angles,ratio,1);
km=b(1);
yy = polyval(b,angles);
figure
plot(angles,ratio,angles,yy)
legend('鉴角曲线','拟合鉴角曲线','Location','NorthEast')
xlabel('\theta_{\it{t}} / degree')
ylabel('\Delta/\Sigma')

%% 单个时刻的信号测角
load Signal
signal = reshape(Signal(:,1,:),[Nrn M]);
data = signal(2401,:);                      % 取一个时刻的信号

% 和波束扫描，寻找相应最大值点作为预计目标角度，即等信号轴
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
disp('单个时刻的测角结果')
angle_fin = angle0 - angle_t

%% 采用5个时刻的信号进行测角,
result = zeros(1,9);
for k = 1:9
    data = signal(2399+k,:);                      % 取一个时刻的信号

    % 和波束扫描，寻找相应最大值点作为预计目标角度，即等信号轴
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
disp('多个时刻平均的测角结果')
mean(result)
