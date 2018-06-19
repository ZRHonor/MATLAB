%% Code to Generate Taylor Weights
% Arik D. Brown
% Original Code Author: F. W. Hopwood
function[wgt] = Taylor(points,sll,nbar)
r = 10^(abs(sll)/20);%主副瓣幅度比
a = log(r+(r*r-1)^0.5) / pi;
sigma2 = nbar^2/(a*a+(nbar-0.5)^2);%波束展宽因子的平方
%--Compute Fm, the Fourier coefficients of the weight set
for m=1:(nbar-1)
    for n=1:(nbar-1)
        f(n,1)=1-m*m/sigma2/(a*a+(n-0.5)*(n-0.5));
        if n ~= m
            f(n,2)=1/(1-m*m/n/n);
        end
        if n==m
            f(n,2)=1;
        end
    end
    g(1,1)=f(1,1);
    g(1,2)=f(1,2);
    for n=2:(nbar-1)
        g(n,1)=g(n-1,1)*f(n,1);
        g(n,2)=g(n-1,2)*f(n,2);
    end
    F(m)=((-1)^(m+1))/2*g(nbar-1,1)*g(nbar-1,2);
end
jj = [1:points]';
xx = (jj-1+0.5)/points - 1/2; %-- column vector
W = ones(size(jj)); %-- column vector
mm = [1:nbar-1]; %-- row vector
W = W + 2*cos(2*pi*xx*mm)*F';
WPK = 1 + 2*sum(F);
wgt = W / WPK;