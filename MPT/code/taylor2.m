clc;
clear all;

c0=3e8;                            %����  
freq=6e9;                          %Ƶ��     
lambda=c0/freq;                    %�����ɿռ䲨��    
D=10*lambda;                       %���߿ھ�ֱ��  
%L=lz*lambda;                       %�״��ֳ���  
%W=wz*lambda;                       %����ֱ��  
SLL=35;                           %�����ƽ    
e=10^(SLL/20);                    %e������Ϊ���������ƽ��  
A=acosh(e)/pi;                        %A =acosh(e)/pi���븱���e����ϵ�����ӡ�  
np=round(A)+6;                     %�ȵ�ƽ�������,np���ڵ���4*A^2/3+1/2  
k0=2*pi/lambda;                    %���ɿռ䲨��  
%���ǵõ�����ͼ��������Ǹ����ƽ�������ƽ��ֵΪ35dB�������ֵû����0������Ϊû�й�һ��������ͼ:  
%[attachment=35678]  
T=0:0.001:pi/2;             %  ������  
c0=3e8;                            %����  
freq=6e9;                          %Ƶ��  
lambda=c0/freq;                    %�����ɿռ䲨��  
D=10*lambda;                       %���߿ھ�ֱ��  
SLL=35;                           %�����ƽ  
e=10^(SLL/20);                    %e������Ϊ���������ƽ��  
A=acosh(e)/pi;                        %A =acosh(e)/pi���븱���e����ϵ�����ӡ�  
np=round(A)+6;                     %�ȵ�ƽ�������,np���ڵ���4*A^2/3+1/2  
k0=2*pi/lambda;                    %���ɿռ䲨��  
%R=(L^2+W^2/4)/W;  
F=zeros(1,np);  
maxs = np;                         % ����Ҫ����һ�ױ�����������������Ŀ  
mju=zeros(1,np);                   % �����������ĸ��������  
incr = 1.0;                        % ��ⲽ������  
v=1;                               % ����һ�ױ����������ĽǱ�       
%h = v+1.9*v^(1/3)+1;               % һ�ױ�����������һ���Ĳ²�ֵ  
mju(1)=0;                          % һ�ױ�����������������  
                             %mju(2)= fzero(@(x)besselj(v,x*pi),h); %һ�ױ����������ĵ�һ�����  
% һ�ױ����������ĵ�1������������  
for s=2:maxs               
    mju(s) = fzero(@(x)besselj(v,x*pi),mju(s-1)+incr);  
end  
mjunp=fzero(@(x)besselj(v,x*pi),mju(np)+incr);%һ�ױ�����������np�����  
sigma=mjunp/(sqrt(A^2+(np-1/2)^2));%sigmaΪ������ϵ��,mjunpΪ��np  
Fm_numerator=ones(1,np-1);                 %Fm�ķ��Ӳ���  
z=zeros(1,np-1);                          %������������  
%���Fm����  
for m=1:np-1  
    for n=1:np-1  
        z(n)=sigma*(A^2+(n-1/2)^2)^(1/2);  
        Fm_numerator(m)=Fm_numerator(m)*(1-(mju(m+1)/z(n))^2);  
    end  
end  
Fm_denominator=ones(1,np-1);%Fm�ķ�ĸ����  
%���Fm��ĸ  
for m=1:np-1  
    for n=1:np-1  
        if n~= m  
            Fm_denominator(m)=Fm_denominator(m)*(1-(mju(m+1)/mju(n+1))^2);  
        end  
    end  
end  
%���ϵ��Fm��������������  
F(1)=1;%ʵ����F0=1  
for i=2:np  
    F(i)=-besselj(0,mju(i)*pi)*Fm_numerator(i-1)/Fm_denominator(i-1);  
    RRR(i-1)=-besselj(0,mju(i)*pi);  
end  
syms u theta;%p;  
%pΪ�ھ��뾶�ѣ�0�ܦѡ�D/2;gΪ�ѵĺ���g���ѣ���ע��Բ�ھ���Ƕ��޹�g=F(1)*besselj(0,mju(1)*pi*p*2/D)/(besselj(0,mju(1)*pi))^2;  
%uΪԲ�ھ�ͨ�ýǱ�����u=2��a*sin��/�ˣ�  
a=D/2;  
u=2*pi*a*sin(theta)/lambda;  
p=0:D/2/(length(T)-1):D/2;%�ھ��뾶�ѱ���  
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