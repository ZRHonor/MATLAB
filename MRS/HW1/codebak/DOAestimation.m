function [doa_CBF, angle]= DOAestimation(X, M, N, P, R)
    R=X*X'/N;
    angle = -90:0.01:90;
    for i =1:length(angle)
        a = exp(-j*2*pi*0.5*[0:M-1]'*sin(pi*angle(i)/180));
        y_CBF(i) = sqrt(abs(a'*R*a));
%         w = (inv(R)*a)/(a'*inv(R)*a);
%         y_MVDR(i) = sqrt(abs(w'*R*a));
    end
    doa_CBF = ESA(angle, y_CBF, P);
%     doa_MVDR = ESA(angle, y_MVDR, P);
end