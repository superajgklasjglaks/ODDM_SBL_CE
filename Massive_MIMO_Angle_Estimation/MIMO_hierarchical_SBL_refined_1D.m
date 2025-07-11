function [h_hat, aod_hat, aoa_hat] = MIMO_hierarchical_SBL_refined_1D(x_p, x_kp, x_lp, Nr, Nt, N_T, M_T, y_T, refined_aod, refined_aoa, on_flag)
% Massive MIMO 分层1D SBL 细化函数
% 模仿ODDM中的hierarchical_SBL_refined_1D
% 输入:
%   x_p: 导频功率
%   x_kp, x_lp: 导频位置参数
%   Nr, Nt: 接收和发射天线数
%   N_T, M_T: 截断矩阵维度
%   y_T: 接收信号向量
%   refined_aod, refined_aoa: 细化的角度网格
%   on_flag: 网格标志
% 输出:
%   h_hat: 估计的信道系数
%   aod_hat, aoa_hat: 估计的角度

virtual_size = length(refined_aod);
rho = 1e-2;
c = 1e-4;
d = 1e-4;

%% 构建细化网格的测量矩阵
delta_aod = zeros(virtual_size, 1);
delta_aoa = zeros(virtual_size, 1);

% 直接使用细化网格作为基础网格
aod_bar = refined_aod;
aoa_bar = refined_aoa;

% 生成测量矩阵
phi_t = zeros(M_T * N_T, virtual_size);
for i = 1:virtual_size
    for nnn = 0:(N_T-1)
        for mmm = 1:M_T
            matrix_idx = nnn * M_T + mmm;
            phi_t(matrix_idx, i) = x_p * ...
                MIMO_Array_Response_Tx(Nt, nnn, x_kp-1, aod_bar(i)) * ...
                MIMO_Array_Response_Rx(Nr, mmm-1, x_lp-1, aoa_bar(i));
        end
    end
end

% 简化的SBL算法
ksi = 1e-3;
Tmax = 200;  % 减少迭代次数

sigma2_bar = sum(abs(y_T).^2) / (100 * N_T * M_T);
alpha = abs(phi_t.' * y_T);
beta0 = 1 / sigma2_bar;

% 迭代过程
for t = 1:Tmax
    alpha_prev = alpha;
    delta_a = diag(alpha);
    
    % 更新mu和Sigma
    C = (1 / beta0) * eye(M_T * N_T) + phi_t * delta_a * phi_t';
    C_inv = inv(C);
    sigma_hbar = delta_a - delta_a * phi_t' * C_inv * phi_t * delta_a;
    mu_hbar = beta0 * sigma_hbar * phi_t' * y_T;
    
    % 计算gamma1
    gamma1 = 1 - real(diag(sigma_hbar)) ./ (alpha + 1e-16);
    
    % 更新alpha
    alpha_temp = abs(mu_hbar).^2 + real(diag(sigma_hbar));
    alpha = -0.5/rho + sqrt(0.25 / rho^2 + alpha_temp / rho);
    
    % 更新beta0
    resid = y_T - phi_t * mu_hbar;
    res2 = norm(resid, 'fro')^2;
    A_beta0 = res2 + 1 / beta0 * sum(gamma1);
    beta0 = (c - 1 + M_T * N_T) / (d + A_beta0);
    
    % 终止条件
    tolerance = norm(alpha - alpha_prev) / norm(alpha_prev);
    if tolerance <= ksi
        break;
    end
end

%% 输出
h_hat = mu_hbar;
aod_hat = aod_bar;
aoa_hat = aoa_bar;

end