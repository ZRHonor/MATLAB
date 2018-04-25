%16QAM中频调制
%升余弦窗成形，滚降系数0.35，符号率1MSybol/s,中频频率21MHz
%给出各级滤波器的系数，时域/频域响应以及信号经过各滤波器的时域/频域图

%*************************产生QAM基带信号*************************%
%产生伪随机序列PN

N=500; %二进制数据长度，长度不足会造成星座点缺失
x=randint(1,N,2); 

%数据分组,串并变换
x1=x(1:2);
x2=x(3:4); %完成第一组转换
for i=1:(N/4-1)  %完成所有点的转换  
    x1=[x1 (x(i*4+1:i*4+2))];
    x2=[x2 (x(i*4+3:i*4+4))];
end

%二-十进制转换（00-0,01-1,10-2,11-3)
xi=x1(1)*2+x1(2);
xq=x2(1)*2+x2(2); %完成第一组转换
n=length(x1);

for i=1:n/2-1    %完成所有点的转换  
    xi=[xi (x1(i*2+1)*2+x1(i*2+2))];
    xq=[xq (x2(i*2+1)*2+x2(i*2+2))];
end

%信号映射
for i=1:n/2
    switch(xi(i))
        case 0
            xi(i)=-3;
        case 1
            xi(i)=-1;
        case 2
            xi(i)=1;
        case 3
            xi(i)=3;
    end
    
    switch(xq(i))
        case 0
            xq(i)=-3;
        case 1
            xq(i)=-1;
        case 2
            xq(i)=1;
        case 3 
            xq(i)=3;
    end
end

%绘制基带信号图
figure(1);
plot(xi,xq,'*'),title('16AQM基带信号星座图');
axis([-4,4,-4,4]);


%***************************脉冲成形*************************************%

%升余弦窗，符号率1MSym/s，设采样频率Fs=4MHz
Fs=4000000; %采样频率设为符号率4倍
t=1/Fs;     %周期
R=0.35;     %滚降系数0.35
N_T=7;      %窗长，设置为标量7
RATE=4;     %采样率
p=rcosfir(R,N_T,RATE,t,'sqrt'); %Design a raised cosine FIR filter，设计升余弦窗

figure(2);  %绘制升余弦窗时域图
stem(p),title('升余弦窗时域图');
figure(3);
freqz(p),title('升余弦窗频域图');
 
%I路和Q路信号成形
xi=upsample(xi,4); %信号内插
xq=upsample(xq,4);
yi=conv(p,xi);  %I路信号与升余弦窗卷积
yq=conv(p,xq);  %Q路信号与升余弦窗卷积
figure(4);      %绘制成形信号图
subplot(2,1,1);
plot(yi),title('I路成形信号');
subplot(2,1,2);
plot(yq),title('Q路成形信号');
figure(5);
plot(db(freqz(yi)/max(abs(freqz(yi))))),title('I路成形信号频谱'),xlabel('频率/Hz'),ylabel('幅度/db');
figure(6);
plot(db(freqz(yq)/max(abs(freqz(yq))))),title('Q路成形信号频谱'),xlabel('频率/Hz'),ylabel('幅度/db');

%********************************HBF********************************%
%设计半带滤波器，阶数为15。
B=firhalfband(14,0.01,'dev'); %15阶半带滤波器时域响应
figure(7);
stem(B);
title('半带滤波器时域图');
figure(8);
plot(db(abs(freqz(B))));
title('半带滤波器频域图');


%I、Q两路信号分别与HBF卷积
%yi=upsample(yi,4);
%yq=upsample(yq,4);
yih=conv(B,yi);
yqh=conv(B,yq);
figure(9);   %绘制经过HBF后的信号
subplot(2,1,1);
plot(yih),title('HBF滤波后I路信号'); %时域信号
subplot(2,1,2);
plot(yqh),title('HBF滤波后Q路信号');
figure(10);
%freqz(yih);
plot(db(freqz(yih)/max(abs(freqz(yih))))),title('HBF滤波后I路信号频谱'),xlabel('频率/Hz'),ylabel('幅度/db') %频域信号
figure(11);
plot(db(freqz(yqh)/max(abs(freqz(yqh))))),title('HBF滤波后Q路信号频谱'),xlabel('频率/Hz'),ylabel('幅度/db');

%*****************************CIC******************************%
%设计CIC
Rc=3; %interpolation factor,内插系数
Mc=1; %defferential delay,延时
N=5;  %number of sections 单元数量
Hm=mfilt.cicinterp(Rc,Mc,N); 
%figure(12);
%stem(Hm);
%I、Q两路信号分别经过CIC滤波
y_fi=filter(Hm,yih);
y_fq=filter(Hm,yqh);
yih=double(yih); yic=double(y_fi); yic=yic/max(abs(yic));%数据类型转换，并归一化
yqh=double(yqh); yqc=double(y_fq); yqc=yqc/max(abs(yqc));
figure(13);
subplot(2,1,1);
plot(yic),title('I路信号经过CIC后的时域响应'),axis([1,900,-1.5,1.5]);
subplot(2,1,2);
plot(yqc),title('Q路信号经过CIC后的时域响应'),axis([1,900,-1.5,1.5]);
figure(14);
%freqz(yic);
plot(db(freqz(yic)/max(abs(freqz(yic))))),title('I路信号经过CIC后的频域响应'),xlabel('频率/Hz'),ylabel('幅度/db');
figure(15);
plot(db(freqz(yqc)/max(abs(freqz(yqc))))),title('Q路信号经过CIC后的频域响应'),xlabel('频率/Hz'),ylabel('幅度/db');


%****************************中频信号*****************************%
%中频21MHz
fc=21000000;
f0=2000000*32;
l=1:length(yic);   %信号序列长度
%F=cos(2*pi*l*fc/f0);  %载波信号
%figure(16);
%plot(abs(freqz(F)));
%figure(17);
%plot(F);
Yi=yic.*cos(2*pi*l*fc/f0);
Yq=yqc.*sin(2*pi*l*fc/f0);
figure(18);
subplot(2,1,1);
plot(Yi/max(abs(Yi))),title('I路中频信号时域图');
subplot(2,1,2);
plot(Yq/max(abs(Yq))),title('q路中频信号时域图');
figure(19)
plot(db(fft(Yi)/max(abs(fft(Yi))))),title('I路中频信号频域图');
figure(20);
plot(db(fft(Yq)/max(abs(fft(Yq))))),title('Q路中频信号频域图');
%I、Q两路信号合成
Y=Yi+Yq;
figure(21);
plot(Y),title('中频输出信号');
figure(23);
plot(abs(fft(Y))),title('中频信号频域图');


