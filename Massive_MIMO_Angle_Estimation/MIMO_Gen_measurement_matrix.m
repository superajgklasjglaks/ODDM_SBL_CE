function [phi_trunc, phi_t, phi_t_aod, phi_t_aoa] = MIMO_Gen_measurement_matrix(x_p, x_kp, x_lp, Nr, Nt, M_T, N_T, M_aoa, N_aod, aod_bar, aoa_bar, delta_aod, delta_aoa)
% Massive MIMO 测量矩阵生成函数
% 模仿ODDM中的Gen_measurement_matrix
% 输入:
%   x_p: 导频功率
%   x_kp, x_lp: 导频位置参数
%   Nr, Nt: 接收和发射天线数
%   M_T, N_T: 截断矩阵维度
%   M_aoa, N_aod: 角度网格维度
%   aod_bar, aoa_bar: 角度网格
%   delta_aod, delta_aoa: 角度偏移
% 输出:
%   phi_trunc: 截断测量矩阵
%   phi_t: 基础测量矩阵
%   phi_t_aod: 发射角导数矩阵
%   phi_t_aoa: 接收角导数矩阵

phi_t = zeros(M_T * N_T, M_aoa * N_aod);
phi_t_aod = zeros(M_T * N_T, M_aoa * N_aod);
phi_t_aoa = zeros(M_T * N_T, M_aoa * N_aod);

for nn = 0:(N_aod-1)
    for mm = 1:M_aoa
        for nnn = 0:(N_T-1)
            for mmm = 1:M_T
                grid_idx = nn * M_aoa + mm;
                matrix_idx = nnn * M_T + mmm;
                
                % 基础测量矩阵
                phi_t(matrix_idx, grid_idx) = x_p * ...
                    MIMO_Array_Response_Tx(Nt, nnn, x_kp-1, aod_bar(grid_idx)) * ...
                    MIMO_Array_Response_Rx(Nr, mmm-1, x_lp-1, aoa_bar(grid_idx));
                
                % 发射角导数矩阵
                phi_t_aod(matrix_idx, grid_idx) = x_p * ...
                    MIMO_wt_derivation_aod(Nt, nnn, x_kp-1, aod_bar(grid_idx)) * ...
                    MIMO_Array_Response_Rx(Nr, mmm-1, x_lp-1, aoa_bar(grid_idx));
                
                % 接收角导数矩阵
                phi_t_aoa(matrix_idx, grid_idx) = x_p * ...
                    MIMO_Array_Response_Tx(Nt, nnn, x_kp-1, aod_bar(grid_idx)) * ...
                    MIMO_wt_derivation_aoa(Nr, mmm-1, x_lp-1, aoa_bar(grid_idx));
            end
        end
    end
end

phi_trunc = phi_t + phi_t_aod * diag(delta_aod) + phi_t_aoa * diag(delta_aoa);

end