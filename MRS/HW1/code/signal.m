function X = signal(M, DOA, N, SNR, QAM)
% 建立信源模型，得到接收的信号
% M 阵元个数
% DOA 信源方向
% N 采样点数
% SNR 信噪比
% QAM QAM调制位数
    P = length(DOA);
    A=exp(-j*2*pi*0.5*[0:M-1]'*sin(DOA));
    %信源模型建立
    for k=1:P
        symbol = randi([0, QAM-1], 1, N);
        S(k,:) = qammod(symbol, QAM);
    end
    X = awgn(A*S,SNR,'measured');

end