%=======================================
% 生成GK101 任意波数据文件的程序
% Copyright GINGKO @2014
% www.cnblogs.com/xiaomagee/p/3930321.html
%=======================================

function arb(x,fre,amp,offs)
%清理工作台
%close all
%clc
%建立文件
fid=fopen('arb0.wvt','wt') ; 
%写入注释头
fprintf(fid,';==== arb file generator for gk101 ====\n');
%写入频率信息
fprintf(fid,'[frequency] = %f;\n',fre);
%写入幅度信息
fprintf(fid,'[amplitude] = %f;\n',amp);
%写入直流偏置信息
fprintf(fid,'[offset] = %f;\n',offs);
%写入日期信息
fprintf(fid,'[date] = %d-%d-%d;\n',year(now),month(now),day(now));
%提取矩阵大、小范围，并写入
fprintf(fid,'[datarange] = %f,%f;\n',min(x),max(x));
[m,n]=size(x);
%提取矩阵长度，并写入
fprintf(fid,'[length] = %d;\n',n);
%写入数据
fprintf(fid,'[data] = \n');
fprintf(fid,'%f,%f,%f,%f,%f,\n',x);     %输出矩阵
%关闭文件
fclose(fid)