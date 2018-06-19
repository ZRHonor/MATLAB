%initialize all the data needed
SNR = 10;
M = 10;						%the number of the array units
d = 0.5;					%the distance of the units normalized with lamda
sigma = 1/sqrt(SNR);		%the varience of the gassian noise
N = 4000;					%the number of the sample times
theta1 = 5;
theta2 = 45;
theta3 = 65;
M_all=(4:1:12);

%calculate the SINR with the M changing	
%%%%%%%%%%%%%%%%% Matched Filter %%%%%%%%%%%%%%%%%%%%%%%
%calculate the spatial steering vectors
As_theta1 = exp(i*2*pi*0.5*sin(theta1*pi/180)*[0:M-1]');
As_theta2 = exp(i*2*pi*0.5*sin(theta2*pi/180)*[0:M-1]');
As_theta3 = exp(i*2*pi*0.5*sin(theta3*pi/180)*[0:M-1]');
%the coherence matrix
S2=eye(M,M);
W_MF=S2^(-1)*conj(As_theta2)/((As_theta2).'*(S2^(-1)*conj(As_theta2)));
    
%%%%%%%%%%%%%%%%% Matched Filter %%%%%%%%%%%%%%%%%%%%%%%R
%tcalculate he coherence matrix
S1=conj(As_theta1)*As_theta1.';
S3=conj(As_theta3)*As_theta3.';
n=normrnd(0,sigma,M,N);
Cov_n=zeros(M,M);
for m=1:N
    Cov_n=Cov_n+n(:,m)*n(:,m)';
end
Cov_n = Cov_n/N;
S_I = Cov_n+(S1+S3);
W_MVDR=S_I^(-1)*conj(As_theta2)/((As_theta2).'*(S_I^(-1)*conj(As_theta2)));
   

SINR=10*log10((W_MF'*conj(As_theta2)*As_theta2.'*W_MF)/(W_MF'*S_I*W_MF))
SINR_MVDR=10*log10((W_MVDR'*conj(As_theta2)*As_theta2.'*W_MVDR)/(W_MVDR'*S_I*W_MVDR))