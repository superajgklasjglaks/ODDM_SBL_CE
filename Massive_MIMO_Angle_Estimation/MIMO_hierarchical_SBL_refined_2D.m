function [H_opt, aod_opt, aoa_opt] = MIMO_hierarchical_SBL_refined_2D(x_p, x_kp, x_lp_factor, x_lp, Nr, Nt, N_T, M_T, y_trunc, aod_refined, aoa_refined, r_aod_fine, r_aoa_fine, aod_max, aoa_max, on_flag)
% Massive MIMO 分层2D SBL 细化函数
% 模仿ODDM中的hierarchical_SBL_refined_2D
% 输入:
%   x_p: 导频功率
%   x_kp, x_lp_factor, x_lp: 导频位置参数
%   Nr, Nt: 接收和发射天线数
%   N_T, M_T: 截断矩阵维度
%   y_trunc: 截断接收信号矩阵
%   aod_refined, aoa_refined: 细化的角度网格
%   r_aod_fine, r_aoa_fine: 细网格分辨率
%   aod_max, aoa_max: 角度最大值
%   on_flag: 网格标志
% 输出:
%   H_opt: 估计的2D信道矩阵
%   aod_opt, aoa_opt: 估计的角度向量

% 将y_trunc重新整形为向量
y_T = reshape(y_trunc.', N_T * M_T, 1);

% 调用分层1D SBL作为基础
[h_hat_temp, aod_hat_temp, aoa_hat_temp] = MIMO_hierarchical_SBL_refined_1D(x_p, x_kp, x_lp, Nr, Nt, N_T, M_T, y_T, aod_refined, aoa_refined, on_flag);

% 计算网格维度
N_aod_refined = length(aod_refined);
M_aoa_refined = length(aoa_refined);

% 将结果重新整形为2D格式
if N_aod_refined * M_aoa_refined == length(h_hat_temp)
    H_opt = reshape(h_hat_temp, M_aoa_refined, N_aod_refined);
else
    % 如果维度不匹配，创建一个合适大小的矩阵
    H_opt = zeros(M_aoa_refined, N_aod_refined);
    % 将非零系数放在合适的位置
    [~, max_idx] = max(abs(h_hat_temp));
    if ~isempty(max_idx)
        [row_idx, col_idx] = ind2sub([M_aoa_refined, N_aod_refined], max_idx);
        H_opt(row_idx, col_idx) = h_hat_temp(max_idx);
    end
end

aod_opt = aod_hat_temp;
aoa_opt = aoa_hat_temp;

end