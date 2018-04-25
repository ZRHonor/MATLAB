function [doa, angle, y]= DOAestimation(X, M, N, P)
    R=X*X'/N;
    angle = -90:0.01:90;
    for i =1:length(angle)
        a = exp(-j*2*pi*0.5*[0:M-1]'*sin(pi*angle(i)/180));
        y(i) = sqrt(abs(a'*R*a));
    end
    doa = ESA(angle, y, P);
end