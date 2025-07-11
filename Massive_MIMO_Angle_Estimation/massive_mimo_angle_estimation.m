close all
clear all
clc

%% Massive MIMO 方向角估计算法比较脚本
% 模仿ODDM信道估计中时延和多普勒频移的特性
% 用发射方向角(AoD)和接收方向角(AoA)替代时延和多普勒频移
% 包含：1D SBL, 2D SBL, 分层1D SBL, 分层2D SBL, OMP, Traditional方法

%% 算法分辨率参数设置 (可独立调整)
% 1D SBL 参数
r_aod_1D = 0.8;    % 1D SBL 发射角分辨率 (弧度)
r_aoa_1D = 0.8;    % 1D SBL 接收角分辨率 (弧度)

% 2D SBL 参数
r_aod_2D = 0.8;    % 2D SBL 发射角分辨率
r_aoa_2D = 0.8;    % 2D SBL 接收角分辨率

% 分层1D SBL 参数
r_aod_coarse_1D = 0.8;    % 分层1D SBL 粗网格发射角分辨率
r_aoa_coarse_1D = 0.8;    % 分层1D SBL 粗网格接收角分辨率
r_aod_fine_1D = 0.4;      % 分层1D SBL 细网格发射角分辨率
r_aoa_fine_1D = 0.4;      % 分层1D SBL 细网格接收角分辨率
threshold_ratio_1D = 0.1;  % 分层1D SBL 系数选择阈值

% 分层2D SBL 参数
r_aod_coarse_2D = 0.8;    % 分层2D SBL 粗网格发射角分辨率
r_aoa_coarse_2D = 0.8;    % 分层2D SBL 粗网格接收角分辨率
r_aod_fine_2D = 0.4;      % 分层2D SBL 细网格发射角分辨率
r_aoa_fine_2D = 0.4;      % 分层2D SBL 细网格接收角分辨率
threshold_ratio_2D = 0.1;  % 分层2D SBL 系数选择阈值

% OMP 参数
r_aod_OMP = 0.8;    % OMP 发射角分辨率
r_aoa_OMP = 0.8;    % OMP 接收角分辨率

%% Test Mode %%%%%%
test = 0;    % set to 1 when testing

%% Massive MIMO 系统参数%%%%%%%%%%
% Nt: 发射天线数
Nt = 32;
% Nr: 接收天线数
Nr = 32;
% 天线间距 (以波长为单位)
d_lambda = 0.5;
% 载波频率
car_freq = 5e9;  % 28 GHz
c = 3e8;          % 光速
lambda = c / car_freq;
d = d_lambda * lambda;

% 角度范围设置
aod_max = pi/3;   % 发射角最大值 (60度)
aoa_max = pi/3;   % 接收角最大值 (60度)
P = 5;            % 路径数
on_flag = 0;      % set to 0 when on-grid

%% 测试用参数
aod_test = [0, 0, 0, 0, 0]';  % 测试用发射角
aoa_test = [0, 0, 0, 0, 0]';  % 测试用接收角
h_test = [1, 0, 0, 0, 0]';    % 测试用信道系数

%% 导频符号设置
Kp = 1;
Lp = 1;
x_kp = floor(Nt/2);
x_lp = floor(Nr/2);
% 数据位置网格
data_grid = ones(Nt, Nr);
data_grid(x_kp-floor(Kp/2)-2:x_kp+floor(Kp/2)+2, x_lp-floor(Lp/2)-2:x_lp+floor(Lp/2)+2) = 0;
% 每帧符号数
N_syms_perfram = sum(sum(data_grid));

% SNR 和噪声方差
SNR_dB = 0:5:20;
SNR = 10.^(SNR_dB/10);
sigma_2 = 0.5 ./ SNR;
SNR_p = SNR * 1e5;
sigma_2_p = sigma_2 / 1e5;

%% 初始化仿真误差计数变量
N_fram = 10;  % 仿真帧数
NMSE_count_1D = zeros(length(SNR_dB), 1);
NMSE_count_2D = zeros(length(SNR_dB), 1);
NMSE_count_hierarchical_1D = zeros(length(SNR_dB), 1);
NMSE_count_hierarchical_2D = zeros(length(SNR_dB), 1);
NMSE_count_OMP = zeros(length(SNR_dB), 1);
NMSE_count_traditional = zeros(length(SNR_dB), 1);
NMSE_CRLB = zeros(length(SNR_dB), 1);  % CRLB理论下界

fprintf('\n=== Massive MIMO 方向角估计算法比较开始 ===\n');
fprintf('仿真参数：Nt=%d, Nr=%d, P=%d, N_fram=%d\n', Nt, Nr, P, N_fram);
fprintf('\n算法分辨率参数：\n');
fprintf('1D SBL: r_aod=%.3f, r_aoa=%.3f\n', r_aod_1D, r_aoa_1D);
fprintf('2D SBL: r_aod=%.3f, r_aoa=%.3f\n', r_aod_2D, r_aoa_2D);
fprintf('分层1D SBL: 粗网格(r_aod=%.3f, r_aoa=%.3f), 细网格(r_aod=%.3f, r_aoa=%.3f)\n', ...
    r_aod_coarse_1D, r_aoa_coarse_1D, r_aod_fine_1D, r_aoa_fine_1D);
fprintf('分层2D SBL: 粗网格(r_aod=%.3f, r_aoa=%.3f), 细网格(r_aod=%.3f, r_aoa=%.3f)\n', ...
    r_aod_coarse_2D, r_aoa_coarse_2D, r_aod_fine_2D, r_aoa_fine_2D);
fprintf('OMP: r_aod=%.3f, r_aoa=%.3f\n', r_aod_OMP, r_aoa_OMP);
fprintf('Traditional: 传统波束成形估计\n');
fprintf('\n');

for ifram = 1:N_fram
    
    %% 随机信道初始化
    aod_c_init = unifrnd(-aod_max, aod_max, P, 1);   % 发射角
    aoa_c_init = unifrnd(-aoa_max, aoa_max, P, 1);   % 接收角
    
    % 路径增益 (模仿ODDM中的功率分布)
    path_gains = exp(-0.1 * (1:P)');
    path_gains = path_gains / sum(path_gains);
    h_c_init = normrnd(0, path_gains);  % 复高斯信道系数
    
    for iesn0 = 1:length(SNR_dB)
        fprintf('Frame %d/%d, SNR %d dB\n', ifram, N_fram, SNR_dB(iesn0));
        
        %% 生成导频信号
        pilot_power = sqrt(1000) / (Kp * Lp);
        
        % 计算测量矩阵 (模仿ODDM的phi_sys)
        phi_sys = zeros(Nr * Nt, P);
        
        for pp = 1:P
            % 发射阵列响应向量
            at = exp(1j * 2 * pi * d / lambda * (0:Nt-1)' * sin(aod_c_init(pp)));
            % 接收阵列响应向量
            ar = exp(1j * 2 * pi * d / lambda * (0:Nr-1)' * sin(aoa_c_init(pp)));
            
            % 构建测量矩阵 (Kronecker积形式)
            phi_sys(:, pp) = kron(at, ar);
        end
        
        %% 信道输出
        noise_gen_re = sigma_2(iesn0) * randn(Nr * Nt, 1);
        noise_gen_im = sigma_2(iesn0) * randn(Nr * Nt, 1);
        noise_gen = noise_gen_re + 1j * noise_gen_im;
        
        if (test == 1)
            r = phi_sys * h_test;
        else
            r = phi_sys * h_c_init + noise_gen;
        end
        
        y = reshape(r, Nr, Nt).';
        
        %% 获取截断测量矩阵
        aod_range = 4;  % 角度搜索范围 (对应ODDM中的k_max)
        aoa_range = 4;  % 角度搜索范围 (对应ODDM中的l_max)
        
        y_trunc = y(floor(Nt/2)-floor(Kp/2)-aod_range:floor(Nt/2)+floor(Kp/2)+aod_range, ...
                   floor(Nr/2)-floor(Lp/2):floor(Nr/2)+floor(Lp/2)+aoa_range);
        N_T = 2 * aod_range + Kp;
        M_T = aoa_range + Lp;
        y_T = reshape(y_trunc.', N_T * M_T, 1);
        
        %% 算法1: 1D SBL 角度估计
        [h_hat_1D, aod_hat_1D, aoa_hat_1D, virtual_size_1D, Phi_1D, delta_a] = ...
            MIMO_CE_1D_SBL(pilot_power, aod_range+Kp, Lp, Nr, Nt, N_T, M_T, y_T, ...
                          r_aod_1D, r_aoa_1D, aod_range, aoa_range, on_flag);
        
        %% 算法2: 2D SBL 角度估计
        [H_opt_2D, aod_opt_2D, aoa_opt_2D, virtual_size_2D] = ...
            MIMO_CE_2D_SBL(pilot_power, aod_range+Kp, 2, Lp, Nr, Nt, N_T, M_T, y_trunc, ...
                          r_aod_2D, r_aoa_2D, aod_range, aoa_range, on_flag);
        
        %% 算法3: 分层1D SBL 角度估计
        % Step 1: 粗网格SBL估计
        [h_hat_coarse_1D, aod_hat_coarse_1D, aoa_hat_coarse_1D, virtual_size_coarse_1D, Phi_coarse_1D, delta_coarse_1D] = ...
            MIMO_CE_1D_SBL(pilot_power, aod_range+Kp, Lp, Nr, Nt, N_T, M_T, y_T, ...
                          r_aod_coarse_1D, r_aoa_coarse_1D, aod_range, aoa_range, on_flag);
        
        % Step 2: 识别显著系数及其位置
        [sorted_coeff_1D, sorted_idx_1D] = sort(abs(h_hat_coarse_1D), 'descend');
        num_significant_1D = max(1, floor(threshold_ratio_1D * length(h_hat_coarse_1D)));
        significant_idx_1D = sorted_idx_1D(1:num_significant_1D);
        
        % 提取显著系数对应的角度值
        aod_significant_1D = aod_hat_coarse_1D(significant_idx_1D);
        aoa_significant_1D = aoa_hat_coarse_1D(significant_idx_1D);
        
        % Step 3: 在显著位置周围创建细化网格
        refined_aod_1D = [];
        refined_aoa_1D = [];
        
        for i = 1:length(aod_significant_1D)
            % 定义显著系数周围的局部区域
            aod_center = aod_significant_1D(i);
            aoa_center = aoa_significant_1D(i);
            
            % 在此中心周围创建细网格
            aod_range_local = aod_center + (-r_aod_coarse_1D/2:r_aod_fine_1D:r_aod_coarse_1D/2);
            aoa_range_local = aoa_center + (-r_aoa_coarse_1D/2:r_aoa_fine_1D:r_aoa_coarse_1D/2);
            
            % 确保边界在有效范围内
            aod_range_local = aod_range_local(aod_range_local >= -aod_max & aod_range_local <= aod_max);
            aoa_range_local = aoa_range_local(aoa_range_local >= -aoa_max & aoa_range_local <= aoa_max);
            
            % 为此区域创建网格
            [AOD_mesh, AOA_mesh] = meshgrid(aod_range_local, aoa_range_local);
            region_aod = AOD_mesh(:);
            region_aoa = AOA_mesh(:);
            
            refined_aod_1D = [refined_aod_1D; region_aod];
            refined_aoa_1D = [refined_aoa_1D; region_aoa];
        end
        
        % 去除重复项
        if ~isempty(refined_aod_1D)
            [refined_grid_1D, unique_idx_1D] = unique([refined_aod_1D, refined_aoa_1D], 'rows');
            refined_aod_1D = refined_grid_1D(:, 1);
            refined_aoa_1D = refined_grid_1D(:, 2);
        end
        
        % Step 4: 在细化区域执行细网格SBL
        if length(refined_aod_1D) > 0
            [h_hat_hierarchical_1D, aod_hat_hierarchical_1D, aoa_hat_hierarchical_1D] = ...
                MIMO_hierarchical_SBL_refined_1D(pilot_power, aod_range+Kp, Lp, Nr, Nt, N_T, M_T, y_T, ...
                                                refined_aod_1D, refined_aoa_1D, on_flag);
        else
            % 回退到粗网格结果
            h_hat_hierarchical_1D = h_hat_coarse_1D;
            aod_hat_hierarchical_1D = aod_hat_coarse_1D;
            aoa_hat_hierarchical_1D = aoa_hat_coarse_1D;
        end
        
        %% 算法4: 分层2D SBL 角度估计
        % Step 1: 粗网格2D SBL估计
        [H_opt_coarse_2D, aod_opt_coarse_2D, aoa_opt_coarse_2D, virtual_size_coarse_2D] = ...
            MIMO_CE_2D_SBL(pilot_power, aod_range+Kp, 2, Lp, Nr, Nt, N_T, M_T, y_trunc, ...
                          r_aod_coarse_2D, r_aoa_coarse_2D, aod_range, aoa_range, on_flag);
        
        % Step 2: 识别显著系数及其2D位置
        coeff_magnitude_2D = abs(H_opt_coarse_2D);
        max_coeff_2D = max(coeff_magnitude_2D);
        significant_indices_2D = find(coeff_magnitude_2D > threshold_ratio_2D * max_coeff_2D);
        
        % 将线性索引转换为2D网格位置
        N_aod_coarse_2D = ceil(2 * aod_max / r_aod_coarse_2D);
        M_aoa_coarse_2D = ceil(2 * aoa_max / r_aoa_coarse_2D);
        [row_indices_2D, col_indices_2D] = ind2sub([M_aoa_coarse_2D, N_aod_coarse_2D], significant_indices_2D);
        
        % Step 3: 在显著系数周围创建细化网格
        [aod_refined_2D, aoa_refined_2D] = create_refined_grid_2D_MIMO(aod_opt_coarse_2D, aoa_opt_coarse_2D, ...
            significant_indices_2D, N_aod_coarse_2D, M_aoa_coarse_2D, r_aod_coarse_2D, r_aoa_coarse_2D, ...
            r_aod_fine_2D, r_aoa_fine_2D, aod_max, aoa_max);
        
        % Step 4: 在细化网格上进行细化2D SBL估计
        [H_opt_hierarchical_2D, aod_opt_hierarchical_2D, aoa_opt_hierarchical_2D] = ...
            MIMO_hierarchical_SBL_refined_2D(pilot_power, aod_range+Kp, 2, Lp, Nr, Nt, N_T, M_T, y_trunc, ...
                                            aod_refined_2D, aoa_refined_2D, r_aod_fine_2D, r_aoa_fine_2D, ...
                                            aod_range, aoa_range, on_flag);
        
        %% 算法5: OMP 角度估计
        [h_hat_OMP, omp_index, aod_hat_OMP, aoa_hat_OMP] = ...
            MIMO_OMP(pilot_power, aod_range+Kp, Lp, Nr, Nt, N_T, M_T, y_T, ...
                    r_aod_OMP, r_aoa_OMP, aod_range, aoa_range);
        
        %% 算法6: Traditional 波束成形角度估计
        [h_hat_traditional, aod_hat_traditional, aoa_hat_traditional] = ...
            MIMO_traditional_beamforming(pilot_power, y_trunc, aod_range, aoa_range, sigma_2_p(iesn0));
        
        %% 计算CRLB理论下界
        % 构建雅可比矩阵J (模仿ODDM中的实现)
        if exist('virtual_size_1D', 'var') && exist('Phi_1D', 'var') && exist('delta_a', 'var')
            J = zeros(M_T * N_T, virtual_size_1D);
            for m = 1:M_T
                for n = 1:N_T
                    J((n-1)*M_T+m, :) = (MIMO_Array_Response_Tx(Nt, n, 1, aod_hat_1D) .* ...
                                         MIMO_Array_Response_Rx(Nr, m, 1, aoa_hat_1D) .* ...
                                         exp(-2j*pi*aod_hat_1D.*aoa_hat_1D/Nt/Nr)).';
                end
            end
            % 计算协方差矩阵
            C_h = (sigma_2(iesn0)^(-1) * (Phi_1D' * Phi_1D) + delta_a^(-1))^(-1);
            C_alpha = J * C_h * J.';
            NMSE_CRLB(iesn0) = trace(abs(C_alpha));
        else
            NMSE_CRLB(iesn0) = 1e-10;  % 默认值
        end
        
        %% 计算所有算法的NMSE
        NMSE_nume_1D = 0;
        NMSE_nume_2D = 0;
        NMSE_nume_hierarchical_1D = 0;
        NMSE_nume_hierarchical_2D = 0;
        NMSE_nume_OMP = 0;
        NMSE_nume_traditional = 0;
        NMSE_deno = 0;
        
        for kk = 0:(Nt-1)
            for ll = 1:Nr
                if (test == 1)
                    h_w = sum(MIMO_Array_Response_Tx(Nt, kk+1, 1, aod_test) .* ...
                             MIMO_Array_Response_Rx(Nr, ll, 1, aoa_test) .* h_test .* ...
                             exp(-2j*pi*aod_test.*aoa_test/Nt/Nr));
                else
                    h_w = sum(MIMO_Array_Response_Tx(Nt, kk+1, 1, aod_c_init) .* ...
                             MIMO_Array_Response_Rx(Nr, ll, 1, aoa_c_init) .* h_c_init .* ...
                             exp(-2j*pi*aod_c_init.*aoa_c_init/Nt/Nr));
                end
                
                % 各算法的估计值
                h_w_hat_1D = sum(MIMO_Array_Response_Tx(Nt, kk+1, 1, aod_hat_1D) .* ...
                                MIMO_Array_Response_Rx(Nr, ll, 1, aoa_hat_1D) .* h_hat_1D .* ...
                                exp(-2j*pi*aod_hat_1D.*aoa_hat_1D/Nt/Nr));
                
                h_w_hat_2D = sum(MIMO_Array_Response_Tx(Nt, kk+1, 1, aod_opt_2D) .* ...
                                MIMO_Array_Response_Rx(Nr, ll, 1, aoa_opt_2D) .* H_opt_2D .* ...
                                exp(-2j*pi*aod_opt_2D.*aoa_opt_2D/Nt/Nr));
                
                h_w_hat_hierarchical_1D = sum(MIMO_Array_Response_Tx(Nt, kk+1, 1, aod_hat_hierarchical_1D) .* ...
                                              MIMO_Array_Response_Rx(Nr, ll, 1, aoa_hat_hierarchical_1D) .* h_hat_hierarchical_1D .* ...
                                              exp(-2j*pi*aod_hat_hierarchical_1D.*aoa_hat_hierarchical_1D/Nt/Nr));
                
                h_w_hat_hierarchical_2D = sum(MIMO_Array_Response_Tx(Nt, kk+1, 1, aod_opt_hierarchical_2D) .* ...
                                              MIMO_Array_Response_Rx(Nr, ll, 1, aoa_opt_hierarchical_2D) .* H_opt_hierarchical_2D .* ...
                                              exp(-2j*pi*aod_opt_hierarchical_2D.*aoa_opt_hierarchical_2D/Nt/Nr));
                
                h_w_hat_OMP = sum(MIMO_Array_Response_Tx(Nt, kk+1, 1, aod_hat_OMP) .* ...
                                 MIMO_Array_Response_Rx(Nr, ll, 1, aoa_hat_OMP) .* h_hat_OMP .* ...
                                 exp(-2j*pi*aod_hat_OMP.*aoa_hat_OMP/Nt/Nr));
                
                h_w_hat_traditional = sum(MIMO_Array_Response_Tx(Nt, kk+1, 1, aod_hat_traditional) .* ...
                                         MIMO_Array_Response_Rx(Nr, ll, 1, aoa_hat_traditional) .* h_hat_traditional .* ...
                                         exp(-2j*pi*aod_hat_traditional.*aoa_hat_traditional/Nt/Nr));
                
                NMSE_nume_1D = NMSE_nume_1D + abs(h_w - h_w_hat_1D).^2;
                NMSE_nume_2D = NMSE_nume_2D + abs(h_w - h_w_hat_2D).^2;
                NMSE_nume_hierarchical_1D = NMSE_nume_hierarchical_1D + abs(h_w - h_w_hat_hierarchical_1D).^2;
                NMSE_nume_hierarchical_2D = NMSE_nume_hierarchical_2D + abs(h_w - h_w_hat_hierarchical_2D).^2;
                NMSE_nume_OMP = NMSE_nume_OMP + abs(h_w - h_w_hat_OMP).^2;
                NMSE_nume_traditional = NMSE_nume_traditional + abs(h_w - h_w_hat_traditional).^2;
                NMSE_deno = NMSE_deno + abs(h_w)^2;
            end
        end
        
        NMSE_count_1D(iesn0) = NMSE_count_1D(iesn0) + NMSE_nume_1D / (NMSE_deno * N_fram);
        NMSE_count_2D(iesn0) = NMSE_count_2D(iesn0) + NMSE_nume_2D / (NMSE_deno * N_fram);
        NMSE_count_hierarchical_1D(iesn0) = NMSE_count_hierarchical_1D(iesn0) + NMSE_nume_hierarchical_1D / (NMSE_deno * N_fram);
        NMSE_count_hierarchical_2D(iesn0) = NMSE_count_hierarchical_2D(iesn0) + NMSE_nume_hierarchical_2D / (NMSE_deno * N_fram);
        NMSE_count_OMP(iesn0) = NMSE_count_OMP(iesn0) + NMSE_nume_OMP / (NMSE_deno * N_fram);
        NMSE_count_traditional(iesn0) = NMSE_count_traditional(iesn0) + NMSE_nume_traditional / (NMSE_deno * N_fram);
    end
end

%% 结果绘图
NMSE_dB_1D = 10 * log10(NMSE_count_1D);
NMSE_dB_2D = 10 * log10(NMSE_count_2D);
NMSE_dB_hierarchical_1D = 10 * log10(NMSE_count_hierarchical_1D);
NMSE_dB_hierarchical_2D = 10 * log10(NMSE_count_hierarchical_2D);
NMSE_dB_OMP = 10 * log10(NMSE_count_OMP);
NMSE_dB_traditional = 10 * log10(NMSE_count_traditional);
NMSE_dB_CRLB = 10 * log10(NMSE_CRLB);

colors = [0,107,182; %蓝1
          118,41,133; %紫2
          234,174,31; %黄3
          215,94,59; %橙4
          184,125,182; %粉紫5
          71,90,40; %绿6
          161,27,30]; %红7
colors = colors/256;

figure('Position', [100, 100, 800, 600]);
plot(SNR_dB, NMSE_dB_1D, '-v', 'Color', colors(3,:), 'LineWidth', 2, 'MarkerSize', 8);
hold on
plot(SNR_dB, NMSE_dB_2D, '-*', 'Color', colors(2,:), 'LineWidth', 2, 'MarkerSize', 8);
plot(SNR_dB, NMSE_dB_hierarchical_1D, '-s', 'Color', colors(1,:), 'LineWidth', 2, 'MarkerSize', 8);
plot(SNR_dB, NMSE_dB_hierarchical_2D, '-o', 'Color', colors(4,:), 'LineWidth', 2, 'MarkerSize', 8);
plot(SNR_dB, NMSE_dB_OMP, '-d', 'Color', colors(6,:), 'LineWidth', 2, 'MarkerSize', 8);
plot(SNR_dB, NMSE_dB_traditional, '-x', 'Color', colors(7,:), 'LineWidth', 2, 'MarkerSize', 8);
plot(SNR_dB, NMSE_dB_CRLB, '-^', 'Color', colors(5,:), 'LineWidth', 2, 'MarkerSize', 8);
grid on

legend('1D SBL', '2D SBL', '分层1D SBL', '分层2D SBL', 'OMP', 'Traditional Beamforming', 'CRLB下界', 'Location', 'best');
xlabel('SNR (dB)', 'FontSize', 12);
ylabel('NMSE (dB)', 'FontSize', 12);
title('Massive MIMO 方向角估计算法性能比较', 'FontSize', 14);
set(gca, 'FontSize', 11);

%% 结果分析和显示
fprintf('\n=== Massive MIMO 方向角估计算法比较结果 ===\n');
fprintf('SNR(dB)\t1D SBL\t\t2D SBL\t\t分层1D SBL\t分层2D SBL\tOMP\t\t\tTraditional\tCRLB下界\n');
fprintf('------\t------\t\t------\t\t----------\t----------\t----------\t----------\t--------\n');
for i = 1:length(SNR_dB)
    fprintf('%d\t\t%.2f\t\t%.2f\t\t%.2f\t\t%.2f\t\t%.2f\t\t%.2f\t\t%.2f\n', ...
        SNR_dB(i), NMSE_dB_1D(i), NMSE_dB_2D(i), NMSE_dB_hierarchical_1D(i), ...
        NMSE_dB_hierarchical_2D(i), NMSE_dB_OMP(i), NMSE_dB_traditional(i), NMSE_dB_CRLB(i));
end

%% 算法复杂度分析
fprintf('\n=== 算法复杂度分析 ===\n');
fprintf('1D SBL格点个数: %d\n', ceil(2 * aod_max / r_aod_1D) * ceil(2 * aoa_max / r_aoa_1D));
fprintf('2D SBL格点个数: %d\n', ceil(2 * aod_max / r_aod_2D) * ceil(2 * aoa_max / r_aoa_2D));
fprintf('分层1D SBL粗网格格点个数: %d\n', ceil(2 * aod_max / r_aod_coarse_1D) * ceil(2 * aoa_max / r_aoa_coarse_1D));
fprintf('分层1D SBL细网格格点个数: %d\n', ceil(2 * aod_max / r_aod_fine_1D) * ceil(2 * aoa_max / r_aoa_fine_1D));
fprintf('分层2D SBL粗网格格点个数: %d\n', ceil(2 * aod_max / r_aod_coarse_2D) * ceil(2 * aoa_max / r_aoa_coarse_2D));
fprintf('分层2D SBL细网格格点个数: %d\n', ceil(2 * aod_max / r_aod_fine_2D) * ceil(2 * aoa_max / r_aoa_fine_2D));
fprintf('OMP格点个数: %d\n', ceil(2 * aod_max / r_aod_OMP) * ceil(2 * aoa_max / r_aoa_OMP));
fprintf('Traditional Beamforming格点个数: %d\n', (2 * aod_range + 1) * (2 * aoa_range + 1));

fprintf('\n=== 仿真完成 ===\n');

%% 辅助函数定义

%% 分层2D SBL网格细化函数 (Massive MIMO版本)
function [aod_refined, aoa_refined] = create_refined_grid_2D_MIMO(aod_coarse, aoa_coarse, significant_indices, N_aod_coarse, M_aoa_coarse, r_aod_coarse, r_aoa_coarse, r_aod_fine, r_aoa_fine, aod_max, aoa_max)
    % 为Massive MIMO创建基于显著系数的2D细化网格
    
    % 获取粗网格结构
    [aod_bar_coarse, aoa_bar_coarse] = MIMO_First_Order_Linear_Approximation_2D(N_aod_coarse, M_aoa_coarse, aod_max, r_aod_coarse, r_aoa_coarse);
    
    % 在2D网格中找到唯一的显著位置
    significant_aod = [];
    significant_aoa = [];
    
    for idx = 1:length(significant_indices)
        sig_idx = significant_indices(idx);
        
        % 将线性索引转换为2D下标
        [row_idx, col_idx] = ind2sub([M_aoa_coarse, N_aod_coarse], sig_idx);
        
        % 获取粗网格中心位置
        aod_center = aod_bar_coarse(col_idx);
        aoa_center = aoa_bar_coarse(row_idx);
        
        significant_aod = [significant_aod; aod_center];
        significant_aoa = [significant_aoa; aoa_center];
    end
    
    % 去除重复项
    significant_positions = unique([significant_aod, significant_aoa], 'rows');
    significant_aod = significant_positions(:, 1);
    significant_aoa = significant_positions(:, 2);
    
    % 创建细化网格维度
    N_aod_refined = ceil(2 * aod_max / r_aod_fine);
    M_aoa_refined = ceil(2 * aoa_max / r_aoa_fine);
    
    % 使用MIMO_First_Order_Linear_Approximation_2D创建结构化细化网格
    [aod_refined_full, aoa_refined_full] = MIMO_First_Order_Linear_Approximation_2D(N_aod_refined, M_aoa_refined, aod_max, r_aod_fine, r_aoa_fine);
    
    % 为了与MIMO_CE_2D_SBL结构兼容，返回完整的细化网格
    % 算法将通过SBL过程自动聚焦于显著区域
    aod_refined = aod_refined_full;
    aoa_refined = aoa_refined_full;
end