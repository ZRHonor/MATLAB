clc;
clear all;
close all;

%------------------------面阵阵---------2009-6-10----------------------------
c=3e8;                     % 光速
fc=35e9;                    % 工作频率（hz）
lamda=c/fc;                % 波长 wave length

N=13;                      % 阵列数
d=0.5*lamda;                 % 阵元间距
L=N*d;                     % 天线长
k=(2*pi)/lamda;            % 波数
fs=5;                      % 采样频率
theta=-90:1/fs:90-1/fs;    % 方位角度范围（采样范围）(rad)
Ns=length(theta);          % 采样点数

theta_d=-0;                % 主波束方向
                           % theta_d=0时w1=w2=...=wM=1,此时即为静态方向图                           
u=sin(theta*pi/180)-sin(theta_d*pi/180);

%-----------------direct--------------------------------------------------
% v=k*d*u;
% 
% E=zeros(1,Ns);
% for n=1:N
%    E=E+sin((2*n-1)/2*v); 
% end
% 
% S=abs(E)/max(abs(E));
% G=20*log10(S);
% 
% figure
% plot(theta,G,'b');axis([-90 90 -60,0]);hold on
% %-----------------------------------------
% v=L*u/numda;
% E=zeros(1,Ns);
% E=(pi*v).*cos(pi*v)./((v-0.5).*(v+0.5))
% S=abs(E)/max(abs(E));
% G=20*log10(S);
% 
% plot(theta,G,'k:');axis([-90 90 -60,0]);hold on
 %----------------------Bayliss分布-----------------------------------
 SLR=-30
 R=10.^(-SLR/20);
 nn=6
 v=L*u/lamda;

 if SLR==-15
     A=1.0079;
     Omega=[1.5124,2.2561,3.1693,4.1264];
 elseif SLR==-20
     A=1.2247;
     Omega=[1.6962,2.2398,3.2473,4.1854];
 elseif SLR==-25
     A=1.4355;
     Omega=[1.8826,2.4943,3.3351,4.2527];
 elseif SLR==-30
     A=1.6413;
     Omega=[2.0708,2.6275,3.4314,4.3276];
 elseif SLR==-35
     A=1.8413;
     Omega=[2.2602,2.7675,3.5352,4.4093];
 elseif SLR==-40
     A=2.0415;
     Omega=[2.4504,2.9123,3.6452,4.4973];      
 end
 
vn=zeros(1,nn-1);
for n=0:(nn-1)
    if n==0
        vn(n+1)=0;
    elseif (n>=1)&&(n<=4)
        vn(n+1)=(nn+0.5)*sqrt(((Omega(n)).^2)/(A.^2+nn.^2)); 
    else
        vn(n+1)=(nn+0.5)*sqrt((A.^2+n.^2)/(A.^2+nn.^2)); 
    end
end

T1=ones(1,Ns);
T2=ones(1,Ns);
T2=(1-(2*v).^2);
for n=1:(nn-1)
    T1=T1.*(1-(v/(vn(n+1))).^2);
    T2=T2.*(1-(v/(n+0.5)).^2);
end
E=zeros(1,Ns);
E=(pi*v).*cos(pi*v).*T1./T2;

S=abs(E)/max(abs(E));
G=20*log10(S);

plot(theta,G,'k');axis([-90 90 -60,0]);hold on
legend('direct','Bayliss');
%-------------(3.89)--------Bayliss分布----------P82----------------

nn=ceil(2*A.^2+0.5)+1
v=L*u/lamda;
vn=zeros(1,nn-1);
for n=0:(nn-1)
    if n==0
        vn(n+1)=0;
    elseif (n>=1)&&(n<=4)
        vn(n+1)=Omega(n); 
    else
        vn(n+1)=sqrt(A.^2+n.^2); 
    end
end

vnn=sqrt(A.^2+nn.^2);
deta=(nn+0.5)/vnn;
T1=ones(1,Ns);
T2=ones(1,Ns);
T2=(1-(2*v).^2);
for n=1:(nn-1)
    vn(n)=sqrt(A.^2+n.^2);
    T1=T1.*(1-(v/(deta*vn(n+1))).^2);
    T2=T2.*(1-(v/(n+0.5)).^2);
end
E=zeros(1,Ns);
E=(pi*v).*cos(pi*v).*T1./T2

S=abs(E)/max(abs(E));
G=20*log10(S);

plot(theta,G,'c--');axis([-90 90 -60,0]);hold on


