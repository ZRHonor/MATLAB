function I=Contrast(g)
%------��������ľ��󣬼�������ĶԱȶȵ�ֵ
%%
g=abs(g);
[N M]=size(g);
mean_g=mean(mean(g.^2));
std_g=std2(g.^2);
I=std_g/mean_g;