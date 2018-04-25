function dataRc = RangeCompression(data, fr, Kr, Nan)
	H_RangeComp=exp(-j*pi*(fr.^2)/Kr)*ones(1,Nan);
	S_f_tm_RangeComp=fftshift(fft(data,[],1),1).*H_RangeComp;
	s_t_tm_RangeComp=ifft(ifftshift(S_f_tm_RangeComp,1),[],1);
	figure;imagesc(abs(s_t_tm_RangeComp));
	xlabel('��λ��');ylabel('������');
    title('������ѹ��Ľ��');
    dataRc = s_t_tm_RangeComp;
end