% function [In] = mybayliss(R, nbar, M)
R = 30;nbar =6;M=16;
%
% 计算bayliss-n分布权系数
%

% 输入参数
% R — 主副瓣比
% n — 等旁瓣个数
% p — 归一化单元位置坐标

% R和n的参考取值如下

% R = 15;                                 %主副瓣比为15dB
% n = 4;                                 

% R = 20;                                 %主副瓣比为20dB
% n = 4;

% R = 25;                                 %主副瓣比为25dB
% n = 5;

% R = 30;                                 %主副瓣比为30dB
% n = 6;

% R = 35;                                 %主副瓣比为35dB
% n = 7;

% R = 40;                                 %主副瓣比为40dB
% n = 8; 

if R == 15
    A = 1.0079;                           %主副瓣比为15dB
    B(1) = 1.5124;
    B(2) = 2.2561;
    B(3) = 3.1693;
    B(4) = 4.1264;
end

if R == 20
    A = 1.2247;                           %主副瓣比为20dB
    B(1) = 1.6962;
    B(2) = 2.3698;
    B(3) = 3.2473;
    B(4) = 4.1854;
end

if R == 25
    A = 1.4355;                           %主副瓣比为25dB
    B(1) = 1.8826;
    B(2) = 2.4943;
    B(3) = 3.3351;
    B(4) = 4.2527;
end

if R == 30
    A = 1.6413;                           %主副瓣比为30dB
    B(1) = 2.0708;
    B(2) = 2.6275;
    B(3) = 3.4314;
    B(4) = 4.3276;
end

if R == 35
    A = 1.8431;                           %主副瓣比为35dB
    B(1) = 2.2602;
    B(2) = 2.7675;
    B(3) = 3.5252;
    B(4) = 4.4093;
end

if R == 40
    A  = 2.0415;                          %主副瓣比为40dB
    B(1) = 2.4504;
    B(2) = 2.9123;
    B(3) = 3.6452;
    B(4) = 4.4973;
end



Phi = zeros(1,nbar);
Phi(2:5) = B(1:4);
for i = 6:nbar
    Phi(i) = sqrt(A^2+(i-1)^2);
end

sigma = (nbar+0.5)/sqrt(A^2 + nbar^2);

b = zeros(1,nbar);

for m = 0:nbar-1
%     multire = 1;
    T1 = 1;
    T2 = 1;
    for n = 1:nbar-1
        T1 = T1*(1-((m+0.5)/(sigma*Phi(n+1)))^2);
    end
    for n = 0:nbar-1
        if n~=m
            T2 = T2*(1-((m+0.5)/(n+0.5))^2);
        end
    end
    b(m+1) = (1/(2*j))*((-1)^m)*(m-0.5)^2*T1/T2;
end
plot(abs(b))

p = linspace(0,1,M);

for i = 1:length(p)
    sumre = 0;
    for n = 0:nbar-1
        sumre = sumre + b(n+1)*sin((pi*p(n+1))/(n+0.5))
    end
    g(i) = sumre;
end
plot(abs(g))


