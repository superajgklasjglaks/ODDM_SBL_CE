function w = MIMO_wt_derivation_aoa(Nr, l, l_in, aoa)
% Massive MIMO 接收角导数函数
% 模仿ODDM中的wt_derivation
% 输入:
%   Nr: 接收天线数
%   l: 当前天线索引
%   l_in: 参考天线索引
%   aoa: 接收方向角 (弧度)
% 输出:
%   w: 角度导数值

% 天线间距 (以波长为单位)
d_lambda = 0.5;

% 计算角度导数
antenna_index = l - l_in;
w = 1j * 2 * pi * d_lambda * antenna_index * cos(aoa) * ...
    exp(1j * 2 * pi * d_lambda * antenna_index * sin(aoa)) / sqrt(Nr);

end