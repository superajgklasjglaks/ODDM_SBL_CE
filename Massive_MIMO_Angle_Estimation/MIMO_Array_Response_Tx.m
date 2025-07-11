function w = MIMO_Array_Response_Tx(Nt, k, k_in, aod)
% Massive MIMO 发射阵列响应函数
% 模仿ODDM中的Sampling_Function_v
% 输入:
%   Nt: 发射天线数
%   k: 当前天线索引
%   k_in: 参考天线索引
%   aod: 发射方向角 (弧度)
% 输出:
%   w: 阵列响应值

% 天线间距 (以波长为单位)
d_lambda = 0.5;

% 计算阵列响应
% 使用均匀线性阵列模型
antenna_index = k - k_in;
w = exp(1j * 2 * pi * d_lambda * antenna_index * sin(aod)) / sqrt(Nt);

end