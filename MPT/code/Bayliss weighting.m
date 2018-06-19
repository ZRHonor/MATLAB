clc;
clear all;
close all;

%------------------------������---------2009-6-10----------------------------
c=3e8;                     % ����
fc=35e9;                    % ����Ƶ�ʣ�hz��
lamda=c/fc;                % ���� wave length

N=13;                      % ������
d=0.5*lamda;                 % ��Ԫ���
L=N*d;                     % ���߳�
k=(2*pi)/lamda;            % ����
fs=5;                      % ����Ƶ��
theta=-90:1/fs:90-1/fs;    % ��λ�Ƕȷ�Χ��������Χ��(rad)
Ns=length(theta);          % ��������

theta_d=-0;                % ����������
                           % theta_d=0ʱw1=w2=...=wM=1,��ʱ��Ϊ��̬����ͼ                           
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
 %----------------------Bayliss�ֲ�-----------------------------------
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
%-------------(3.89)--------Bayliss�ֲ�----------P82----------------

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


