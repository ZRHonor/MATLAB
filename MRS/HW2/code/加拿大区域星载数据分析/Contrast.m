function I=Contrast(g)
%------根据输入的矩阵，计算输出的对比度的值
%%
g=abs(g);
[N M]=size(g);
mean_g=mean(mean(g.^2));
std_g=std2(g.^2);
I=std_g/mean_g;