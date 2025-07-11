function w = MIMO_wt_derivation_aod(Nt, k, k_in, aod)
% Massive MIMO 发射角导数函数
% 模仿ODDM中的wv_derivation
% 输入:
%   Nt: 发射天线数
%   k: 当前天线索引
%   k_in: 参考天线索引
%   aod: 发射方向角 (弧度)
% 输出:
%   w: 角度导数值

% 天线间距 (以波长为单位)
d_lambda = 0.5;

% 计算角度导数
antenna_index = k - k_in;
w = 1j * 2 * pi * d_lambda * antenna_index * cos(aod) * ...
    exp(1j * 2 * pi * d_lambda * antenna_index * sin(aod)) / sqrt(Nt);

end