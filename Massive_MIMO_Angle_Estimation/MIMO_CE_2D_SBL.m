function [H_opt, aod_opt, aoa_opt, virtual_size] = MIMO_CE_2D_SBL(x_p, x_kp, x_lp_factor, x_lp, Nr, Nt, N_T, M_T, y_trunc, r_aod, r_aoa, aod_max, aoa_max, on_flag)
% Massive MIMO 2D SBL 角度估计函数
% 简化版本，模仿ODDM中的CE_2D_SBL
% 输入:
%   x_p: 导频功率
%   x_kp, x_lp_factor, x_lp: 导频位置参数
%   Nr, Nt: 接收和发射天线数
%   N_T, M_T: 截断矩阵维度
%   y_trunc: 截断接收信号矩阵
%   r_aod, r_aoa: 角度分辨率
%   aod_max, aoa_max: 角度最大值
%   on_flag: 网格标志
% 输出:
%   H_opt: 估计的2D信道矩阵
%   aod_opt, aoa_opt: 估计的角度向量
%   virtual_size: 虚拟网格大小

N_aod = ceil(2 * aod_max / r_aod);
M_aoa = ceil(2 * aoa_max / r_aoa);
virtual_size = N_aod * M_aoa;

% 生成角度网格
[aod_bar, aoa_bar] = MIMO_First_Order_Linear_Approximation_2D(N_aod, M_aoa, aod_max, r_aod, r_aoa);

% 简化的2D SBL实现
% 将y_trunc重新整形为向量
y_T = reshape(y_trunc.', N_T * M_T, 1);

% 调用1D SBL作为基础
[h_hat_temp, aod_hat_temp, aoa_hat_temp, ~, ~, ~] = MIMO_CE_1D_SBL(x_p, x_kp, x_lp, Nr, Nt, N_T, M_T, y_T, r_aod, r_aoa, aod_max, aoa_max, on_flag);

% 将结果重新整形为2D格式
H_opt = reshape(h_hat_temp, M_aoa, N_aod);
aod_opt = aod_hat_temp;
aoa_opt = aoa_hat_temp;

end