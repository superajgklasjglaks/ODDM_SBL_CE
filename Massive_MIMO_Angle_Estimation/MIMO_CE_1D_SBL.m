function [h_hat, aod_hat, aoa_hat, virtual_size, phi_trunc, delta_a] = MIMO_CE_1D_SBL(x_p, x_kp, x_lp, Nr, Nt, N_T, M_T, y_T, r_aod, r_aoa, aod_max, aoa_max, on_flag)
% Massive MIMO 1D SBL 角度估计函数
% 模仿ODDM中的CE_1D_SBL
% 输入:
%   x_p: 导频功率
%   x_kp, x_lp: 导频位置参数
%   Nr, Nt: 接收和发射天线数
%   N_T, M_T: 截断矩阵维度
%   y_T: 接收信号向量
%   r_aod, r_aoa: 角度分辨率
%   aod_max, aoa_max: 角度最大值
%   on_flag: 网格标志
% 输出:
%   h_hat: 估计的信道系数
%   aod_hat, aoa_hat: 估计的角度
%   virtual_size: 虚拟网格大小
%   phi_trunc: 截断测量矩阵
%   delta_a: 对角矩阵

N_aod = ceil(2 * aod_max / r_aod);   % 发射角虚拟采样网格
M_aoa = ceil(2 * aoa_max / r_aoa);   % 接收角虚拟采样网格
virtual_size = N_aod * M_aoa;
rho = 1e-2;
c = 1e-4;
d = 1e-4;

%% 开始离网格角度估计
delta_aod = zeros(virtual_size, 1);
delta_aoa = zeros(virtual_size, 1);
[aod_bar, aoa_bar] = MIMO_First_Order_Linear_Approximation(N_aod, M_aoa, aod_max, r_aod, r_aoa);
[phi_trunc, phi_t, phi_t_aod, phi_t_aoa] = MIMO_Gen_measurement_matrix(x_p, x_kp, x_lp, Nr, Nt, M_T, N_T, M_aoa, N_aod, aod_bar, aoa_bar, delta_aod, delta_aoa);

ksi = 1e-3;
Tmax = 5e2;

sigma2_bar = sum(abs(y_T).^2) / (100 * N_T * M_T);
P_hat = floor(M_T * N_T / log(virtual_size));
alpha = abs(phi_trunc.' * y_T);
beta0 = 1 / sigma2_bar;

% 迭代过程
for t = 1:Tmax
    phi_trunc = phi_t + phi_t_aod * diag(delta_aod) + phi_t_aoa * diag(delta_aoa);
    alpha_prev = alpha;
    delta_a = diag(alpha);
    
    % 算法步骤3: 更新mu和Sigma
    C = (1 / beta0) * eye(M_T * N_T) + phi_trunc * delta_a * phi_trunc';
    C_inv = inv(C);
    sigma_hbar = delta_a - delta_a * phi_trunc' * C_inv * phi_trunc * delta_a;  % 方程41
    mu_hbar = beta0 * sigma_hbar * phi_trunc' * y_T;    % 方程40
    
    % 算法步骤4: 计算gamma1用于更新beta0
    gamma1 = 1 - real(diag(sigma_hbar)) ./ (alpha + 1e-16);
    
    % 更新alpha
    alpha_temp = abs(mu_hbar).^2 + real(diag(sigma_hbar));
    alpha = -0.5/rho + sqrt(0.25 / rho^2 + alpha_temp / rho);  % 方程49
    
    % 更新beta0
    resid = y_T - phi_trunc * mu_hbar;
    res2 = norm(resid, 'fro')^2;
    A_beta0 = res2 + 1 / beta0 * sum(gamma1); % 方程51
    beta0 = (c - 1 + M_T * N_T) / (d + A_beta0);
    
    % 终止条件
    tolerance = norm(alpha - alpha_prev) / norm(alpha_prev);
    if tolerance <= ksi
        break;
    end
    
    [~, idx] = sort(alpha, 'descend');
    delta_aod_prev = delta_aod;
    delta_aoa_prev = delta_aoa;
    mmph = (mu_hbar * mu_hbar' + sigma_hbar);
    
    if (tolerance < 1000 * ksi)
        %% 更新delta_aod (发射角偏移)
        idx = idx(1:P_hat);
        
        pHpv = phi_t_aod' * phi_t_aod;
        tptlv = (phi_t + phi_t_aoa * diag(delta_aoa_prev));
        A_aod = real(pHpv .* (conj(mu_hbar) * mu_hbar.' + sigma_hbar.')); % 方程54
        b_aod1 = diag(mu_hbar) * phi_t_aod.' * conj(y_T);
        b_aod2 = diag(mmph * tptlv' * phi_t_aod);
        b_aod = real(b_aod1 - b_aod2);    % 方程55
        
        A_T_aod = A_aod(idx, idx);
        b_T_aod = b_aod(idx);
        
        temp = delta_aod;
        delta_T_aod = A_T_aod \ b_T_aod;
        delta_aod_upgrate = zeros(virtual_size, 1);
        delta_aod_upgrate(idx) = temp(idx);
        flag = ((max(svd(A_T_aod)) / min(svd(A_T_aod))) < 1e6);
        
        if (flag == 0)
            for pp = 1:P_hat
                temp_delta_aod = delta_aod(idx);
                temp_delta_aod(pp) = 0;
                if A_T_aod(pp, pp) == 0
                    delta_aod_upgrate(idx(pp)) = 0;
                    continue;
                else
                    delta_aod_upgrate(idx(pp)) = (b_T_aod(pp) - A_T_aod(pp, :) * temp_delta_aod) / A_T_aod(pp, pp);   % 方程57
                end
            end
        else
            delta_aod_upgrate = zeros(virtual_size, 1);
            delta_aod_upgrate(idx) = delta_T_aod;
        end
        
        delta_aod_upgrate(delta_aod_upgrate < (-r_aod)/2) = (-r_aod)/2;
        delta_aod_upgrate(delta_aod_upgrate > (r_aod)/2) = (r_aod)/2;
        if (on_flag == 1)
            delta_aod = delta_aod_upgrate;
        end
        
        %% 更新delta_aoa (接收角偏移)
        pHpt = phi_t_aoa' * phi_t_aoa;
        tptlt = (phi_t + phi_t_aod * diag(delta_aod_prev));
        A_aoa = real(pHpt .* (conj(mu_hbar) * mu_hbar.' + sigma_hbar.')); % 方程54
        b_aoa1 = diag(mu_hbar) * phi_t_aoa.' * conj(y_T);
        b_aoa2 = diag((mmph * tptlt') * phi_t_aoa);
        b_aoa = real(b_aoa1 - b_aoa2);    % 方程55
        
        A_T_aoa = A_aoa(idx, idx);
        b_T_aoa = b_aoa(idx);
        
        temp = delta_aoa;
        delta_T_aoa = A_T_aoa \ b_T_aoa;
        delta_aoa_upgrate = zeros(virtual_size, 1);
        delta_aoa_upgrate(idx) = temp(idx);
        flag = ((max(svd(A_T_aoa)) / min(svd(A_T_aoa))) < 1e6);
        
        if (flag == 0)
            for pp = 1:P_hat
                temp_delta_aoa = delta_aoa(idx);
                temp_delta_aoa(pp) = 0;
                if A_T_aoa(pp, pp) == 0
                    delta_aoa_upgrate(idx(pp)) = 0;
                    continue;
                else
                    delta_aoa_upgrate(idx(pp)) = (b_T_aoa(pp) - A_T_aoa(pp, :) * temp_delta_aoa) / A_T_aoa(pp, pp);   % 方程57
                end
            end
        else
            delta_aoa_upgrate = zeros(virtual_size, 1);
            delta_aoa_upgrate(idx) = delta_T_aoa;
        end
        
        delta_aoa_upgrate(delta_aoa_upgrate < (-r_aoa)/2) = (-r_aoa)/2;
        delta_aoa_upgrate(delta_aoa_upgrate > (r_aoa)/2) = (r_aoa)/2;
        if (on_flag == 1)
            delta_aoa = delta_aoa_upgrate;
        end
    end
end

%% 输出
h_hat = mu_hbar;
aod_hat = aod_bar + delta_aod;
aoa_hat = aoa_bar + delta_aoa;

end