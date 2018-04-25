function I=Entropy(g)
%------根据输入的矩阵，计算输出的熵的值
%%
[N M]=size(g);
S=sum(sum(abs(g).^2));
temp=abs(g).^2/S.*log10(S./abs(g).^2);
I=sum(sum(temp));
