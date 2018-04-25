clear all; close all; clc;

c = 3e8;
fc = 1e9;
lambda = c/fc;

%% 目标区域参数
% 目标区域方位向范围[Xmin,Xmax]
Xmin = 0;
Xmax = 50;

 %成像区域中线
Yc=10000;
% 目标区域距离向范围[Yc-Y0,Yc+Y0]
% 成像宽度为2*Y0
Y0=500;

%% 轨道参数
% SAR的运动速度
V=100;
% 高度
H = 5000;
% 最短距离；
R0 = sqrt(Y_c^2 + H^2)

%% 天线参数
% 方位向天线长度
D=4;
% SAR合成孔径长度

