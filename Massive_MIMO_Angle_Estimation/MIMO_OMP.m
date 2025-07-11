function [h_hat, omp_index, aod_hat, aoa_hat] = MIMO_OMP(x_p, x_kp, x_lp, Nr, Nt, N_T, M_T, y_T, r_aod, r_aoa, aod_max, aoa_max)
% Massive MIMO OMP 角度估计函数
% 模仿ODDM中的OMP
% 输入:
%   x_p: 导频功率
%   x_kp, x_lp: 导频位置参数
%   Nr, Nt: 接收和发射天线数
%   N_T, M_T: 截断矩阵维度
%   y_T: 接收信号向量
%   r_aod, r_aoa: 角度分辨率
%   aod_max, aoa_max: 角度最大值
% 输出:
%   h_hat: 估计的信道系数
%   omp_index: OMP选择的索引
%   aod_hat, aoa_hat: 估计的角度

N_aod = ceil(2 * aod_max / r_aod);
M_aoa = ceil(2 * aoa_max / r_aoa);
virtual_size = N_aod * M_aoa;

% 生成角度网格
[aod_bar, aoa_bar] = MIMO_First_Order_Linear_Approximation(N_aod, M_aoa, aod_max, r_aod, r_aoa);

% 构建字典矩阵
A = zeros(M_T * N_T, virtual_size);
for i = 1:virtual_size
    for nnn = 0:(N_T-1)
        for mmm = 1:M_T
            matrix_idx = nnn * M_T + mmm;
            A(matrix_idx, i) = x_p * ...
                MIMO_Array_Response_Tx(Nt, nnn, x_kp-1, aod_bar(i)) * ...
                MIMO_Array_Response_Rx(Nr, mmm-1, x_lp-1, aoa_bar(i));
        end
    end
end

% OMP算法参数
max_iter = min(10, floor(M_T * N_T / 4));  % 最大迭代次数
threshold = 1e-6;  % 残差阈值

% 初始化
residual = y_T;
selected_indices = [];
selected_atoms = [];

% OMP主循环
for iter = 1:max_iter
    % 计算内积
    correlations = abs(A' * residual);
    
    % 找到最大相关性的原子
    [~, max_idx] = max(correlations);
    
    % 添加到选择集合
    selected_indices = [selected_indices, max_idx];
    selected_atoms = [selected_atoms, A(:, max_idx)];
    
    % 最小二乘求解
    if size(selected_atoms, 2) == 1
        coeffs = selected_atoms \ y_T;
    else
        coeffs = pinv(selected_atoms) * y_T;
    end
    
    % 更新残差
    residual = y_T - selected_atoms * coeffs;
    
    % 检查终止条件
    if norm(residual) < threshold
        break;
    end
end

% 构建稀疏解
h_hat = zeros(virtual_size, 1);
if ~isempty(selected_indices)
    h_hat(selected_indices) = coeffs;
end

% 输出
omp_index = selected_indices;
aod_hat = aod_bar;
aoa_hat = aoa_bar;

end