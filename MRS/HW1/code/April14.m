clear all;
clc;
close all;

M = 10;                         % ��Ԫ��Ŀ
DOA = [5 45 60].*pi/180;		% ��������
d = 0.5;						% ��Ԫ����
N = 400;						% ��������
QAM = 16;						% 16QAM����
SNR = 10;

% ��������A
P = length(DOA);
for i = 1:P
	 A(:,i)=exp(-j*2*pi*[0:M-1]'*sin(DOA(i))*0.5);
end

% ��Դģ��
% TODO QAM���ƣ���������������������������������
for j = 1:P
	s = sind(1:N);
	S(j,:) = s;
end

% �����ź�ģ�ͽ���
 X=awgn(A*S, SNR);

% Э�������
R=X*X'/N;

angle = -90:1:90;
for k = 1:length(angle)
    a = exp(-j*2*pi*0.5*[0:M-1]'*sin(angle(k)));
    y(k) = sqrt(abs(a'*R*a));
end
plot(angle, 10*log(y/max(y)))
axis([-90 90 -80 10]);
xlabel('theta/degree');
ylabel('��һ���ռ���/dB');