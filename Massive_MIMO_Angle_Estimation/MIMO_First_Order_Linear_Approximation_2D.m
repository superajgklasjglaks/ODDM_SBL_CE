function [aod_bar, aoa_bar] = MIMO_First_Order_Linear_Approximation_2D(N_aod, M_aoa, aod_max, r_aod, r_aoa)
% Massive MIMO 2D一阶线性近似函数
% 模仿ODDM中的First_Order_Linear_Approximation_2D
% 输入:
%   N_aod: 发射角网格数
%   M_aoa: 接收角网格数
%   aod_max: 发射角最大值
%   r_aod: 发射角分辨率
%   r_aoa: 接收角分辨率
% 输出:
%   aod_bar: 发射角网格向量
%   aoa_bar: 接收角网格向量

% 生成发射角网格
aod_bar = linspace(-aod_max, aod_max, N_aod);
% 生成接收角网格 (假设接收角范围与发射角相同)
aoa_bar = linspace(-aod_max, aod_max, M_aoa);

end