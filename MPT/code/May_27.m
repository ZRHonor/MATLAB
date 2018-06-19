c = 3e8;        % 光速
T = 400e-6;     % 发射脉宽
B = 1e6;        % 调频带宽
Rmin = 50e3; Rmax = 250e3;  % 距离门
R = [80e3 200e3 80e3];      % 目标位置
RCS = [1 1 1]               % 目标相对有效反射面

R = R - Rmin;               % 目标位置，相对距离门下限
K = B/T;                    % 调频斜率
Rrec = Rmax-Rmin;           % 接收门宽 m
Trec = 2*Rrec/c;            % 接受门宽 s
Fs = 5*B; Ts = 1/Fs;        % 计算机仿真的采样率和周期
N0 = ceil(T/Ts);            % 发射脉宽对应的采样点数
N = ceil(Trec/Ts);          % 接收窗对应的采样点数



