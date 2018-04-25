function x = StatSigGenerate(M,N,DOA,SNR,SignalMode,d)

ld = length(DOA);
if strcmp(SignalMode,'Independent')
	st = randn(ld,N)+1j*randn(ld,N);
elseif strcmp(SignalMode,'Coherent')
    	st = [];
    	st1 = randn(1,N)+1j*randn(1,N);
	for k = 1:ld
		st = [st;st1];
	end
end
st = st/sqrt(trace(st*st'/N)/ld);
nt = randn(M,N)+1j*randn(M,N);
nt = nt/sqrt(trace(nt*nt'/N)/M);

SNR = ones(1,ld)*SNR;
Amp = diag(10.^(SNR/20));
A = exp(1j*2*pi*[0:M-1]'*sind(DOA)*d);
x = A*Amp*st+nt;
end