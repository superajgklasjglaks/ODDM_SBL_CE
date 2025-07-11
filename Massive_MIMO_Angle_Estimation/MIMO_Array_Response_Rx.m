function w = MIMO_Array_Response_Rx(Nr, l, l_in, aoa)
% Massive MIMO 接收阵列响应函数
% 模仿ODDM中的Sampling_Function_t
% 输入:
%   Nr: 接收天线数
%   l: 当前天线索引
%   l_in: 参考天线索引
%   aoa: 接收方向角 (弧度)
% 输出:
%   w: 阵列响应值

% 天线间距 (以波长为单位)
d_lambda = 0.5;

% 计算阵列响应
% 使用均匀线性阵列模型
antenna_index = l - l_in;
w = exp(1j * 2 * pi * d_lambda * antenna_index * sin(aoa)) / sqrt(Nr);

end