function theta_d = ESA(angle, ss, N)
% 估计信源角度
%   angle和ss为空间谱数据
%   N = 信源个数
    [pks,locs] = findpeaks(ss, angle, 'SortStr', 'descend');
    theta_d = sort(locs(1:N));
end