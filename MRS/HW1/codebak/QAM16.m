clc;
clear
close all;
M=64;
k=log2(M);
n=120000;
samp=1;%过采样率？？？？？？？？？
%snr=0:1:14;


x=randi([0 1],1,n);
stem(x(1:50),'filled');
title('二进制随机比特流');
xlabel('比特序列');
ylabel('信号幅度');

x6=reshape(x,k,length(x)/k);
xsym=bi2de(x6.','left-msb');%2到32进制序列

figure;
stem(xsym(1:50));   %画出相应的16进制信号序列    32
title('16进制随机信号');
xlabel('信号序列');
ylabel('信号幅度');

y=modulate(modem.qammod(M),xsym);
scatterplot(y);

text(real(y)+0.1,imag(y),dec2bin(xsym));%?

axis([-8 8 -8 8]);

snr=15+10*log10(k)-10*log10(samp);
%snrtem=10.^(snr/10);%信噪比？？？？
%pn=1./snrtem;
%sigma=sqrt(pn);
yn=awgn(y,snr,'measured');%加入高白噪声
h=scatterplot(yn,samp,0,'b.');%经过信道后含白噪声的星座图
hold on;
scatterplot(y,1,0,'k+',h);% H must be a valid handle to a figure
    % that was previously generated by SCATTERPLOT.  Default for H is [], which
    %  causes SCATTERPLOT to create a new figure.
title('接收信号星座图');
legend('含噪声接收信号','不含噪声信号');
axis([-8 8 -8 8]);
hold on;
%eyediagram(yn,2);              %眼图
yd=demodulate(modem.qamdemod(M),yn);%解调出来的是16进制信号
z=de2bi(yd,'left-msb');
z=reshape(z.',numel(z),1);
[nuber_of_errors,bit_error_rate]=biterr(x,z)
semilogy(snr,bit_error_rate,'ro');


