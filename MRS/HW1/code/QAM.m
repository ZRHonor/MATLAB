%16QAM��Ƶ����
%�����Ҵ����Σ�����ϵ��0.35��������1MSybol/s,��ƵƵ��21MHz
%���������˲�����ϵ����ʱ��/Ƶ����Ӧ�Լ��źž������˲�����ʱ��/Ƶ��ͼ

%*************************����QAM�����ź�*************************%
%����α�������PN

N=500; %���������ݳ��ȣ����Ȳ�������������ȱʧ
x=randint(1,N,2); 

%���ݷ���,�����任
x1=x(1:2);
x2=x(3:4); %��ɵ�һ��ת��
for i=1:(N/4-1)  %������е��ת��  
    x1=[x1 (x(i*4+1:i*4+2))];
    x2=[x2 (x(i*4+3:i*4+4))];
end

%��-ʮ����ת����00-0,01-1,10-2,11-3)
xi=x1(1)*2+x1(2);
xq=x2(1)*2+x2(2); %��ɵ�һ��ת��
n=length(x1);

for i=1:n/2-1    %������е��ת��  
    xi=[xi (x1(i*2+1)*2+x1(i*2+2))];
    xq=[xq (x2(i*2+1)*2+x2(i*2+2))];
end

%�ź�ӳ��
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

%���ƻ����ź�ͼ
figure(1);
plot(xi,xq,'*'),title('16AQM�����ź�����ͼ');
axis([-4,4,-4,4]);


%***************************�������*************************************%

%�����Ҵ���������1MSym/s�������Ƶ��Fs=4MHz
Fs=4000000; %����Ƶ����Ϊ������4��
t=1/Fs;     %����
R=0.35;     %����ϵ��0.35
N_T=7;      %����������Ϊ����7
RATE=4;     %������
p=rcosfir(R,N_T,RATE,t,'sqrt'); %Design a raised cosine FIR filter����������Ҵ�

figure(2);  %���������Ҵ�ʱ��ͼ
stem(p),title('�����Ҵ�ʱ��ͼ');
figure(3);
freqz(p),title('�����Ҵ�Ƶ��ͼ');
 
%I·��Q·�źų���
xi=upsample(xi,4); %�ź��ڲ�
xq=upsample(xq,4);
yi=conv(p,xi);  %I·�ź��������Ҵ����
yq=conv(p,xq);  %Q·�ź��������Ҵ����
figure(4);      %���Ƴ����ź�ͼ
subplot(2,1,1);
plot(yi),title('I·�����ź�');
subplot(2,1,2);
plot(yq),title('Q·�����ź�');
figure(5);
plot(db(freqz(yi)/max(abs(freqz(yi))))),title('I·�����ź�Ƶ��'),xlabel('Ƶ��/Hz'),ylabel('����/db');
figure(6);
plot(db(freqz(yq)/max(abs(freqz(yq))))),title('Q·�����ź�Ƶ��'),xlabel('Ƶ��/Hz'),ylabel('����/db');

%********************************HBF********************************%
%��ư���˲���������Ϊ15��
B=firhalfband(14,0.01,'dev'); %15�װ���˲���ʱ����Ӧ
figure(7);
stem(B);
title('����˲���ʱ��ͼ');
figure(8);
plot(db(abs(freqz(B))));
title('����˲���Ƶ��ͼ');


%I��Q��·�źŷֱ���HBF���
%yi=upsample(yi,4);
%yq=upsample(yq,4);
yih=conv(B,yi);
yqh=conv(B,yq);
figure(9);   %���ƾ���HBF����ź�
subplot(2,1,1);
plot(yih),title('HBF�˲���I·�ź�'); %ʱ���ź�
subplot(2,1,2);
plot(yqh),title('HBF�˲���Q·�ź�');
figure(10);
%freqz(yih);
plot(db(freqz(yih)/max(abs(freqz(yih))))),title('HBF�˲���I·�ź�Ƶ��'),xlabel('Ƶ��/Hz'),ylabel('����/db') %Ƶ���ź�
figure(11);
plot(db(freqz(yqh)/max(abs(freqz(yqh))))),title('HBF�˲���Q·�ź�Ƶ��'),xlabel('Ƶ��/Hz'),ylabel('����/db');

%*****************************CIC******************************%
%���CIC
Rc=3; %interpolation factor,�ڲ�ϵ��
Mc=1; %defferential delay,��ʱ
N=5;  %number of sections ��Ԫ����
Hm=mfilt.cicinterp(Rc,Mc,N); 
%figure(12);
%stem(Hm);
%I��Q��·�źŷֱ𾭹�CIC�˲�
y_fi=filter(Hm,yih);
y_fq=filter(Hm,yqh);
yih=double(yih); yic=double(y_fi); yic=yic/max(abs(yic));%��������ת��������һ��
yqh=double(yqh); yqc=double(y_fq); yqc=yqc/max(abs(yqc));
figure(13);
subplot(2,1,1);
plot(yic),title('I·�źž���CIC���ʱ����Ӧ'),axis([1,900,-1.5,1.5]);
subplot(2,1,2);
plot(yqc),title('Q·�źž���CIC���ʱ����Ӧ'),axis([1,900,-1.5,1.5]);
figure(14);
%freqz(yic);
plot(db(freqz(yic)/max(abs(freqz(yic))))),title('I·�źž���CIC���Ƶ����Ӧ'),xlabel('Ƶ��/Hz'),ylabel('����/db');
figure(15);
plot(db(freqz(yqc)/max(abs(freqz(yqc))))),title('Q·�źž���CIC���Ƶ����Ӧ'),xlabel('Ƶ��/Hz'),ylabel('����/db');


%****************************��Ƶ�ź�*****************************%
%��Ƶ21MHz
fc=21000000;
f0=2000000*32;
l=1:length(yic);   %�ź����г���
%F=cos(2*pi*l*fc/f0);  %�ز��ź�
%figure(16);
%plot(abs(freqz(F)));
%figure(17);
%plot(F);
Yi=yic.*cos(2*pi*l*fc/f0);
Yq=yqc.*sin(2*pi*l*fc/f0);
figure(18);
subplot(2,1,1);
plot(Yi/max(abs(Yi))),title('I·��Ƶ�ź�ʱ��ͼ');
subplot(2,1,2);
plot(Yq/max(abs(Yq))),title('q·��Ƶ�ź�ʱ��ͼ');
figure(19)
plot(db(fft(Yi)/max(abs(fft(Yi))))),title('I·��Ƶ�ź�Ƶ��ͼ');
figure(20);
plot(db(fft(Yq)/max(abs(fft(Yq))))),title('Q·��Ƶ�ź�Ƶ��ͼ');
%I��Q��·�źźϳ�
Y=Yi+Yq;
figure(21);
plot(Y),title('��Ƶ����ź�');
figure(23);
plot(abs(fft(Y))),title('��Ƶ�ź�Ƶ��ͼ');


