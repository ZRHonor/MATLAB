function [In] = bayliss_n(R, n, p)
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

sigma = (n + 0.5)/sqrt(A^2 + n^2);        %展宽因子

for k = 1:n-1
    if k <= 4
        psi_n(k) = B(k);                  %方向图的零点位置
    else
        psi_n(k) = sqrt(A^2 + k^2);
    end
end

% b1 = zeros(1,length(p));
% for i = 0:n-1
%     multi = 1;
%     for k = 1:n-1
%         if k==i
%             temp = (1-((i+0.5)/(sigma*psi_n(k)))^2);
%             multi = multi*temp;
%         else
%             temp = (1-((i+0.5)/(sigma*psi_n(k)))^2)/(1-((i+0.5)/(k+0.5))^2);
%             multi = multi*temp;
%         end
%     end
%     b1(i+1) = ((i-0.5)^2*multi*(-1)^i)/(2*j);
% end
% b1
% 
% g = zeros(1,length(p));
% n = [1:length(p)]-1;
% 
% for i = 1:length(p)
%     sumrel = 0;
%     for k = 0:n-1
%         temp = b(k+1)*sin((pi*p(i+1))/(k+0.5));
%         sumrel = sumrel +temp;
%     end
%     g(i) = sumrel;
% end
% 
% In = g;
for k = 0:n-1 
    m = 0:n-1;
    i = find(m ~= k);
    m = m(i);    
    b(k+1) = pi*(k+0.5)*cos(pi*(k+0.5))*prod(1-(k+0.5)^2/sigma^2./psi_n.^2)./prod(1-(k+0.5)^2./(m+0.5).^2);
end

% In = ones(1,p);
for k = 1:length(p)
    In(k) = sum(b.*sin(pi*p(k).*((0:n-1)+0.5)));
end

In = -In/max(abs(In));
In = In*(length(p)/sum(abs(In)));
for i=1:floor(length(p/2))
    In(i) = -In(length(p)+1-i);
end
