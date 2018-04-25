%=================================================
%=================================================
function signal=cs(RawSignal,Rc,Rref,r0)                       %reference range,r0?????
%=================================================
Vs=116*2;               %????
c=3.0e+8;              %??
WaveLength=0.03;       %??
fc=0;
fs=160e6;%30.0e+6;
ChirpRate=15.55e6/WaveLength;%5.0e+11;     %????k         这是什么东西
PulseLength=9e-6;%50.0e-6;   %??????Tr
PRF= 750/2;%2100;              %???
N_range=1024*2%512*4;
N_azimuth=1024;
N_signal=ceil((PulseLength)*fs);
%===============Chirp Scaling=====================
%=================================================
%1.Azimuth FFT
S=fft(RawSignal,N_azimuth,1);

%2.Chirp Scaling Phase Multiply
Rref=Rc;
TOU=(-(N_range-1)/2:(N_range-1)/2)/fs;
%TOU is a vector, the first element is -(N_range-1)/(2*fs),the last element is (N_range-1)/(2*fs)
Fs=[(0:N_azimuth/2-1),(-N_azimuth/2:-1)]/N_azimuth*PRF; %Fs is a vector
Cs=1./sqrt(1-(WaveLength*Fs/(2*Vs)).^2)-1;
TOUref=2*Rref*(1.+Cs)/c;
Kstmp=(WaveLength*Fs/(2*Vs)).^2;
Ks=ChirpRate./( 1+ChirpRate*Rref*2*WaveLength/(c^2)*( Kstmp./((1-Kstmp).^(3/2))));
Ft=[(0:N_range/2-1),(-N_range/2:-1)]/N_range*fs;

TOU0=2*r0/c;
for n=1:N_azimuth
    Phase2=-j*pi*Ks(n)*2*Cs(n)*TOU0*TOU+j*pi*Ks(n)*(-TOU0^2*Cs(n)+4*TOU0*Rref*Cs(n)*(1+Cs(n))/c);
    S(n,:)=S(n,:).*exp(-j*pi*Ks(n)*Cs(n)*((TOU-TOUref(n)).^2)+Phase2);
    %S(n,:)=S(n,:).*exp(-j*pi*Ks(n)*Cs(n)*((TOU-TOUref(n)).^2));
end

%3.Range FFT
S=fft(S,N_range,2);

%4.RCMC,Range Compreesion, and SRC Multiply
for n=1:N_azimuth
    S(n,:)=S(n,:).*exp(-j*pi*(Ft.^2)./(Ks(n)*(1+Cs(n)))).*exp(j*4*pi/c*Rref*Ft*Cs(n));
end

%5.Range IFFT
S=ifft(S,N_range,2);

%6.Azimuth Filter and Phase Residual
TOU=2*Rc/c;
for n=1:N_azimuth
    Phas1=-j*2*pi/WaveLength*c*TOU*(1-sqrt(1-(WaveLength*Fs(n)/(2*Vs))^2));
    Phas2=+j*4*pi/c^2*Ks(n)*(1+Cs(n))*Cs(n)*(Rc-Rref)^2;
    S(n,:)=S(n,:).*exp(Phas1+Phas2);
end

%7.Azimuth IFFT
signal=ifft(S,N_azimuth,1);