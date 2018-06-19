clc;
clear all;

c0=3e8;                            %光速  
freq=6e9;                          %频率     
lambda=c0/freq;                    %λ自由空间波长    
D=10*lambda;                       %天线口径直径  
%L=lz*lambda;                       %雷达罩长度  
%W=wz*lambda;                       %底面直径  
SLL=35;                           %副瓣电平    
e=10^(SLL/20);                    %e，即ξ为主、副瓣电平比  
A=acosh(e)/pi;                        %A =acosh(e)/pi是与副瓣比e相联系的因子。  
np=round(A)+6;                     %等电平副瓣个数,np大于等于4*A^2/3+1/2  
k0=2*pi/lambda;                    %自由空间波数  
%可是得到方向图结果并不是副瓣电平与主瓣电平差值为35dB，主瓣峰值没有在0处，因为没有归一化吗？如下图:  
%[attachment=35678]  
T=0:0.001:pi/2;             %  俯仰角  
c0=3e8;                            %光速  
freq=6e9;                          %频率  
lambda=c0/freq;                    %λ自由空间波长  
D=10*lambda;                       %天线口径直径  
SLL=35;                           %副瓣电平  
e=10^(SLL/20);                    %e，即ξ为主、副瓣电平比  
A=acosh(e)/pi;                        %A =acosh(e)/pi是与副瓣比e相联系的因子。  
np=round(A)+6;                     %等电平副瓣个数,np大于等于4*A^2/3+1/2  
k0=2*pi/lambda;                    %自由空间波数  
%R=(L^2+W^2/4)/W;  
F=zeros(1,np);  
maxs = np;                         % 所需要计算一阶贝塞尔函数的零点的数目  
mju=zeros(1,np);                   % 贝塞尔函数的根，即零点  
incr = 1.0;                        % 求解步进常数  
v=1;                               % 代表一阶贝塞尔函数的角标       
%h = v+1.9*v^(1/3)+1;               % 一阶贝塞尔函数第一零点的猜测值  
mju(1)=0;                          % 一阶贝塞尔函数第零个零点  
                             %mju(2)= fzero(@(x)besselj(v,x*pi),h); %一阶贝塞尔函数的第一个零点  
% 一阶贝塞尔函数的第1个及后面的零点  
for s=2:maxs               
    mju(s) = fzero(@(x)besselj(v,x*pi),mju(s-1)+incr);  
end  
mjunp=fzero(@(x)besselj(v,x*pi),mju(np)+incr);%一阶贝塞尔函数第np个零点  
sigma=mjunp/(sqrt(A^2+(np-1/2)^2));%sigma为σ膨胀系数,mjunp为μnp  
Fm_numerator=ones(1,np-1);                 %Fm的分子部分  
z=zeros(1,np-1);                          %见程序函数方程  
%求解Fm分子  
for m=1:np-1  
    for n=1:np-1  
        z(n)=sigma*(A^2+(n-1/2)^2)^(1/2);  
        Fm_numerator(m)=Fm_numerator(m)*(1-(mju(m+1)/z(n))^2);  
    end  
end  
Fm_denominator=ones(1,np-1);%Fm的分母部分  
%求解Fm分母  
for m=1:np-1  
    for n=1:np-1  
        if n~= m  
            Fm_denominator(m)=Fm_denominator(m)*(1-(mju(m+1)/mju(n+1))^2);  
        end  
    end  
end  
%求解系数Fm，见程序函数方程  
F(1)=1;%实际是F0=1  
for i=2:np  
    F(i)=-besselj(0,mju(i)*pi)*Fm_numerator(i-1)/Fm_denominator(i-1);  
    RRR(i-1)=-besselj(0,mju(i)*pi);  
end  
syms u theta;%p;  
%p为口径半径ρ；0≤ρ≤D/2;g为ρ的函数g（ρ）；注：圆口径与角度无关g=F(1)*besselj(0,mju(1)*pi*p*2/D)/(besselj(0,mju(1)*pi))^2;  
%u为圆口径通用角变量，u=2πa*sinθ/λ；  
a=D/2;  
u=2*pi*a*sin(theta)/lambda;  
p=0:D/2/(length(T)-1):D/2;%口径半径ρ变量  
g=zeros(1,length(p));  
gg=zeros(1,np);  
%g(1)=F(1)*besselj(0,mju(1)*pi*p(1)*2/D)/(besselj(0,mju(1)*pi))^2;  
for i=1:length(g)  
    gg(1)=F(1)*besselj(0,mju(1)*pi*p(i)*2/D)/(besselj(0,mju(1)*pi))^2;  
    for m=2:np  
        gg(m)=gg(m-1)+F(m)*besselj(0,mju(m)*pi*p(i)*2/D)/(besselj(0,mju(m)*pi))^2;  
        g(i)=gg(m);  
    end  
    end  
for i=1:length(T)  
    Q(i)=quad(@(p)2.*pi.*g(i).*besselj(0,2.*pi.*freq./c0.*sin(T(i))).*p,0,10*c0/freq/2);  
end  
plot(T,20*log10(abs(Q)));%figure;%hold on;  