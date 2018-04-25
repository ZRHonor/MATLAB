% RD 算法----SAR点目标仿真
% 不考虑 RCM 
% 2006年12月
%-----------------------------------------------
% 参数设定 单位都统一到 m s Hz
clear;clc;
R0=10.e3;       % m 
Vr=100;         % m/s
Tr=10.e-6;      % s
B=50.e6;        % Hz
f0=10.e9;       % Hz
deltafd=50;     % Hz
Fr=100.e6;      % Hz
Fa=100;         % Hz
Naz=128; 
Nrg=2048;
c=3.e8;
lambda=c/f0;    % 波长
La=1;
Kr=B/Tr;        % Hz/s  LFM 的调频斜率
% 距离向取 9.5 km--10.5 km 的范围 
Rrmin=9500;Rrmax=10500;     % m
Ntarget=3;
Rtarget=[10000;10050;10050];Raz=[0;0;30];       % 分别为目标横向和纵向的距离 Raz 为距参考点 fd=0 的距离
% 求与目标有关的参数
tp=-Tr/2:1/Fr:Tr/2;  
Ntr=length(tp);   
ta=linspace(-Naz/Fa/2,Naz/Fa/2,Naz);
% 二维参数
N_st=ceil(2*(Rtarget-Rrmin)/c*Fr)+1;% N_st N_end 分别为点目标回波起始和结束时的采样点
N_end=N_st+Ntr-1;             % 考虑 RCM 时不准确
% 存储回波
S0=zeros(Naz,Nrg);        % 回波信号  
for k=1:Ntarget
    phase=zeros(Naz,Nrg);
    Rta(k,:)=sqrt(Rtarget(k)^2+(Vr*ta-Raz(k)).^2);   % Rta 为斜距矩阵
    A=((tp.^2).'*ones(1,Naz)).';
    B=(Rta(k,:).'*ones(1,Ntr));
    phase(:,N_st(k):N_end(k))= pi*Kr*A-4*pi*f0/c*B;                          
    S0(:,N_st(k):N_end(k))=S0(:,N_st(k):N_end(k))+exp(j*phase(:,N_st(k):N_end(k)));     %sinc(a(k)*tka(k,m))^2*
end  
% 距离压缩
S1=zeros(Naz,Nrg);
S2=S1;
s2=zeros(Naz,Nrg);
S1=(fft(S0.')).';
hrg=exp(-j*pi*Kr*tp.^2);
Hrg=fft(hrg,Nrg);
H=(Hrg.'*ones(1,Naz)).';
S2=S1.*H;
s2=circshift(ifft(S2.').',[0,-(Ntr+1)]); 
%dbs2=20*log10(abs(s2)/max(max(abs(s2))));
% 方位向压缩
Ta=deltafd*lambda*R0/(2*Vr^2);
ta2=1/Fa:1/Fa:Ta;
Ka=2*Vr^2/(lambda*R0);
S3=fft(s2);
haz=exp(j*pi*Ka*ta2.^2);
Haz=fft(haz,Naz).'*ones(1,Nrg);
S4=S3.*Haz;
s4=ifft(S4);
%dbs4=20*log10(abs(s4)/max(max(abs(s4))));
% 绘图
Max=max(max(abs(s4)));
Min=min(min(abs(s4)));
s6=round(255*(abs(s4)-Min)/(Max-Min));
colormap(gray)
image(s6);
