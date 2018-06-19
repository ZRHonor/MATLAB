function WeightVector = getWeightVector(M, sll, nbar, theta0)
% M 阵元个数
% sll 旁瓣幅度dB
% nbar the first nb nulls of S is the same with the dieal space factor ,and n1>=2*A^2+0.5
    A = acosh(10^(abs(sll/20)))/pi;
    nbar = floor(2*A^2+0.5)+2;
    
    n1 = 1:nbar-1;
    sigma = nbar/sqrt(A^2+(nbar-0.5)^2);
    u0 = sigma*sqrt(A^2 + (n1-0.5)^2);
    
    
    
end