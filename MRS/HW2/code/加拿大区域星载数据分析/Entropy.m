function I=Entropy(g)
%------��������ľ��󣬼���������ص�ֵ
%%
[N M]=size(g);
S=sum(sum(abs(g).^2));
temp=abs(g).^2/S.*log10(S./abs(g).^2);
I=sum(sum(temp));
