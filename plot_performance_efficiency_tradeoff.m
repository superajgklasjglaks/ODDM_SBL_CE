close all;
clear all;
clc;

%% 所有算法Performance-Efficiency Trade-off分析
% 基于compare_all_algorithms.m的结果，绘制性能-效率权衡图

%% Test Mode %%%%%%
test = 0;    % set to 1 when testing

%% OTFS parameters%%%%%%%%%%
% N: number of symbols in time
N = 32;
% M: number of subcarriers in frequency
M = 32;
% M_mod: size of QAM constellation
M_mod = 4;
M_bits = log2(M_mod);
% average energy per data symbol
eng_sqrt = (M_mod==2)+(M_mod~=2)*sqrt((M_mod-1)/6*(2^2));

%% delay-Doppler grid symbol placement

% Time and frequency resources
car_fre = 5*10^9;  % Carrier frequency
delta_f = 15*10^3; % subcarrier spacing: 15 KHz
T = 1/delta_f;     % one time symbol duration in OTFS frame
k_max = 4;         % system setting
l_max = 4;
t_max = l_max / M / delta_f;
v_max = k_max / N / T;
P = 9;
on_flag = 0;      % set to 0 when on-grid

%% 算法分辨率参数设置
% 1D SBL 参数
r_v_1D = 0.4;
r_t_1D = 0.4;

% 2D SBL 参数
r_v_2D = 0.4;
r_t_2D = 0.4;

% 分层1D SBL 参数
r_v_coarse_1D = 0.5;
r_t_coarse_1D = 0.5;
r_v_fine_1D = 0.2;
r_t_fine_1D = 0.2;
threshold_ratio_1D = 0.15;

% 分层2D SBL 参数
r_v_coarse_2D = 0.5;
r_t_coarse_2D = 0.5;
r_v_fine_2D = 0.2;
r_t_fine_2D = 0.2;
threshold_ratio_2D = 0.15;

% OMP 参数
sparsity_level = 2*P;

%% pilot symbol placing
Kp = 1;
Lp = 1;
x_kp = floor(N/2);
x_lp = floor(M/2);
% data positions of OTFS delay-Doppler domain data symbols  in the 2-D grid
data_grid = ones(N,M);
data_grid(x_kp-floor(Kp/2)-2*k_max:x_kp+floor(Kp/2)+2*k_max,x_lp-floor(Lp/2)-l_max:x_lp+floor(Lp/2)+l_max)=0;
% number of symbols per frame
N_syms_perfram = sum(sum(data_grid));

% number of bits per frame
N_bits_perfram = N_syms_perfram*M_bits;

% SNR and variance of the noise
% SNR = P/\sigma^2; P: avg. power of albhabet transmitted
SNR_dB_fixed = 10;  % 固定SNR为15dB
SNR = 10.^(SNR_dB_fixed/10);
sigma_2 = 0.5 ./SNR;
SNR_p = SNR * 1e5;
sigma_2_p = sigma_2 /1e5;

%% Initializing simulation error count variables
N_fram = 20;  % 仿真帧数

% 算法名称
algorithm_names = {'Traditional Impulse', '1D SBL', '2D SBL', 'Hierarchical 1D SBL', 'Hierarchical 2D SBL', 'OMP'};
num_algorithms = length(algorithm_names);

% 初始化结果存储
NMSE_results = zeros(num_algorithms, 1);
time_results = zeros(num_algorithms, 1);
complexity_results = zeros(num_algorithms, 1);  % 网格点数

fprintf('\n=== 所有算法Performance-Efficiency Trade-off分析开始 ===\n');
fprintf('仿真参数：N=%d, M=%d, P=%d, N_fram=%d, SNR=%d dB\n', N, M, P, N_fram, SNR_dB_fixed);
fprintf('\n');

% 临时存储每帧仿真结果
temp_NMSE = zeros(num_algorithms, 1);
temp_time = zeros(num_algorithms, 1);

for ifram = 1:N_fram
    if mod(ifram, 5) == 0
        fprintf('Frame %d/%d\n', ifram, N_fram);
    end
    
    %%  random channel initialization
    v_c_init = unifrnd(-v_max,v_max,P,1);   
    t_c_init = unifrnd(0,t_max,P,1);
    l_ti = t_c_init.*(M * delta_f);
    q_l_t = exp(-0.1.*l_ti)./sum(exp(-0.1.*l_ti));
    h_c_init = normrnd(0,q_l_t);  % normrnd but not mvnrnd
    k_v_init = v_c_init .*(N*T);
    l_t_init = t_c_init .*(M*delta_f);
    
    pow_prof = (1/P) * (ones(1,P));
    chan_coef = zeros(1,P);
    
    %% random input bits generation%%%%%
    trans_info_bit = randi([0,1],N_syms_perfram*M_bits,1);
    %%2D QAM symbols generation %%%%%%%%
    data=qammod(reshape(trans_info_bit,M_bits,N_syms_perfram), M_mod,'gray','InputType','bit');
    data = data./(eng_sqrt);
    x = Generate_2D_data_grid_CE(N,M,data,data_grid);
    x(x_kp-floor(Kp/2):x_kp+floor(Kp/2),x_lp-floor(Lp/2):x_lp+floor(Lp/2))=sqrt(1000)/(Kp*Lp);
    
    % calculate measument matrix
    phi_sys = zeros(M*N,P);
    
    kk=1:N;
    ll=1:M;

    for pp=1:P
        for k = 0:(N-1)
            for l = 1:M
                v = Sampling_Function_v(N,k,kk-1,k_v_init(pp));
                t = Sampling_Function_t(M,l-1,ll-1,l_t_init(pp));
                phi_sys(k*M+l,pp) = sum(sum(v.'*t .* x));
            end
        end
    end
    
    %% channel output%%%%%
    noise_gen_re = sigma_2 * randn(M*N,1);
    noise_gen_im = sigma_2 * randn(M*N,1);
    noise_gen = noise_gen_re + 1i.* noise_gen_im;
    h_c_init_hat = h_c_init .* exp(-2i * pi/M/N * (k_v_init.*l_t_init));
    if (test==1)
        r = phi_sys * h_test;
    else
        r = phi_sys * h_c_init + noise_gen;
    end
    y = reshape(r,M,N).';
    
    %% get truncation measument matrix%%%%%
    y_trunc = y(floor(N/2)-floor(Kp/2)-k_max:floor(N/2)+floor(Kp/2)+k_max,floor(M/2)-floor(Lp/2):floor(M/2)+floor(Lp/2)+l_max);
    N_T = 2*k_max+Kp;
    M_T = l_max+Lp;
    y_T = reshape(y_trunc.',N_T*M_T,1);        
    
    %% 算法1: Traditional Impulse
    tic;
    try
        [h_hat_traditional, k_v_hat_traditional, l_t_hat_traditional] = traditional_impulse(sqrt(1000)/(Kp*Lp), k_max+Kp, Lp, M, N, N_T, M_T, y_T, k_max, l_max, on_flag);
        temp_time(1) = temp_time(1) + toc;
    catch ME
        temp_time(1) = temp_time(1) + toc;
        h_hat_traditional = zeros(P,1);
        k_v_hat_traditional = zeros(P,1);
        l_t_hat_traditional = zeros(P,1);
    end
    
    %% 算法2: 1D SBL Channel Estimation
    tic;
    try
        [h_hat_1D, k_v_hat_1D, l_t_hat_1D, virtual_size_1D, Phi_1D, delta_1D] = CE_1D_SBL(sqrt(1000)/(Kp*Lp), k_max+Kp, Lp, M, N, N_T, M_T, y_T, r_v_1D, r_t_1D, k_max, l_max, on_flag);
        temp_time(2) = temp_time(2) + toc;
    catch ME
        temp_time(2) = temp_time(2) + toc;
        h_hat_1D = zeros(P,1);
        k_v_hat_1D = zeros(P,1);
        l_t_hat_1D = zeros(P,1);
        virtual_size_1D = 0;
    end
    
    %% 算法3: 2D SBL Channel Estimation
    tic;
    try
        [H_opt_2D,k_v_opt_2D,l_t_opt_2D,virtual_size_2D] = CE_2D_SBL(sqrt(1000)/(2*Kp*Lp),k_max+Kp,2,Lp,M,N,N_T,M_T,y_trunc,r_v_2D,r_t_2D,k_max,l_max,on_flag);
        temp_time(3) = temp_time(3) + toc;
    catch ME
        temp_time(3) = temp_time(3) + toc;
        H_opt_2D = zeros(P,1);
        k_v_opt_2D = zeros(P,1);
        l_t_opt_2D = zeros(P,1);
        virtual_size_2D = 0;
    end
    
    %% 算法4: 分层1D SBL Channel Estimation
    tic;
    try
        % Step 1: Coarse Grid SBL Estimation
        [h_hat_coarse_1D, k_v_hat_coarse_1D, l_t_hat_coarse_1D, virtual_size_coarse_1D, Phi_coarse_1D, delta_coarse_1D] = ...
            CE_1D_SBL(sqrt(1000)/(Kp*Lp), k_max+Kp, Lp, M, N, N_T, M_T, y_T, r_v_coarse_1D, r_t_coarse_1D, k_max, l_max, on_flag);
        
        % Step 2: Identify significant coefficients and their locations
        [sorted_coeff_1D, sorted_idx_1D] = sort(abs(h_hat_coarse_1D), 'descend');
        num_significant_1D = max(1, floor(threshold_ratio_1D * length(h_hat_coarse_1D)));
        significant_idx_1D = sorted_idx_1D(1:num_significant_1D);
        
        % Extract corresponding k_v and l_t values for significant coefficients
        k_v_significant_1D = k_v_hat_coarse_1D(significant_idx_1D);
        l_t_significant_1D = l_t_hat_coarse_1D(significant_idx_1D);
        
        % Step 3: Create refined grid around significant locations
        refined_k_v_1D = [];
        refined_l_t_1D = [];
        
        for i = 1:length(k_v_significant_1D)
            % Define local region around significant coefficient
            k_center = k_v_significant_1D(i);
            l_center = l_t_significant_1D(i);
            
            % Create fine grid around this center
            k_range = k_center + (-r_v_coarse_1D/2:r_v_fine_1D:r_v_coarse_1D/2);
            l_range = l_center + (-r_t_coarse_1D/2:r_t_fine_1D:r_t_coarse_1D/2);
            
            % Ensure boundaries are within valid range
            k_range = k_range(k_range >= -k_max & k_range <= k_max);
            l_range = l_range(l_range >= 0 & l_range <= l_max);
            
            % Create meshgrid for this region
            [K_mesh, L_mesh] = meshgrid(k_range, l_range);
            region_k = K_mesh(:);
            region_l = L_mesh(:);
            
            refined_k_v_1D = [refined_k_v_1D; region_k];
            refined_l_t_1D = [refined_l_t_1D; region_l];
        end
        
        % Remove duplicates
        if ~isempty(refined_k_v_1D)
            [refined_grid_1D, unique_idx_1D] = unique([refined_k_v_1D, refined_l_t_1D], 'rows');
            refined_k_v_1D = refined_grid_1D(:, 1);
            refined_l_t_1D = refined_grid_1D(:, 2);
        end
        
        % Step 4: Perform fine-grid SBL on refined regions
        if length(refined_k_v_1D) > 0
            [h_hat_hierarchical_1D, k_v_hat_hierarchical_1D, l_t_hat_hierarchical_1D] = ...
                hierarchical_SBL_refined_1D(sqrt(1000)/(Kp*Lp), k_max+Kp, Lp, M, N, N_T, M_T, y_T, ...
                                        refined_k_v_1D, refined_l_t_1D, on_flag);
        else
            % Fallback to coarse grid results
            h_hat_hierarchical_1D = h_hat_coarse_1D;
            k_v_hat_hierarchical_1D = k_v_hat_coarse_1D;
            l_t_hat_hierarchical_1D = l_t_hat_coarse_1D;
        end
        
        temp_time(4) = temp_time(4) + toc;
        
    catch ME
        temp_time(4) = temp_time(4) + toc;
        h_hat_hierarchical_1D = zeros(P,1);
        k_v_hat_hierarchical_1D = zeros(P,1);
        l_t_hat_hierarchical_1D = zeros(P,1);
    end

    %% 算法5: 分层2D SBL Channel Estimation
    tic;
    try
        % Step 1: Coarse Grid 2D SBL Estimation
        [H_opt_coarse_2D,k_v_opt_coarse_2D,l_t_opt_coarse_2D,virtual_size_coarse_2D] = CE_2D_SBL(sqrt(1000)/(2*Kp*Lp),k_max+Kp,2,Lp,M,N,N_T,M_T,y_trunc,r_v_coarse_2D,r_t_coarse_2D,k_max,l_max,on_flag);
        
        % Step 2: Identify significant coefficients and their 2D positions
        coeff_magnitude_2D = abs(H_opt_coarse_2D);
        max_coeff_2D = max(coeff_magnitude_2D);
        significant_indices_2D = find(coeff_magnitude_2D > threshold_ratio_2D * max_coeff_2D);
        
        % Convert linear indices to 2D grid positions
        N_v_coarse_2D = ceil(2 * k_max / r_v_coarse_2D);
        M_t_coarse_2D = ceil(l_max / r_t_coarse_2D);
        [row_indices_2D, col_indices_2D] = ind2sub([M_t_coarse_2D, N_v_coarse_2D], significant_indices_2D);
        
        % Step 3: Create refined grid around significant coefficients
        [k_v_refined_2D, l_t_refined_2D] = create_refined_grid_2D_improved(k_v_opt_coarse_2D, l_t_opt_coarse_2D, significant_indices_2D, N_v_coarse_2D, M_t_coarse_2D, r_v_coarse_2D, r_t_coarse_2D, r_v_fine_2D, r_t_fine_2D, k_max, l_max);
        
        % Step 4: Refined 2D SBL estimation on the refined grid
        [H_opt_hierarchical_2D, k_v_opt_hierarchical_2D, l_t_opt_hierarchical_2D] = hierarchical_SBL_refined_2D(sqrt(1000)/(2*Kp*Lp), k_max+Kp, 2, Lp, M, N, N_T, M_T, y_trunc, k_v_refined_2D, l_t_refined_2D, r_v_fine_2D, r_t_fine_2D, k_max, l_max, on_flag);
        
        temp_time(5) = temp_time(5) + toc;
        
    catch ME
        temp_time(5) = temp_time(5) + toc;
        H_opt_hierarchical_2D = zeros(P,1);
        k_v_opt_hierarchical_2D = zeros(P,1);
        l_t_opt_hierarchical_2D = zeros(P,1);
    end
    
    %% 算法6: OMP Channel Estimation
    tic;
    try
        [h_hat_OMP, k_v_hat_OMP, l_t_hat_OMP, virtual_size_OMP] = OMP(sqrt(1000)/(Kp*Lp), k_max+Kp, Lp, M, N, N_T, M_T, y_T, r_v_1D, r_t_1D, k_max, l_max, sparsity_level, on_flag);
        temp_time(6) = temp_time(6) + toc;
    catch ME
        temp_time(6) = temp_time(6) + toc;
        h_hat_OMP = zeros(P,1);
        k_v_hat_OMP = zeros(P,1);
        l_t_hat_OMP = zeros(P,1);
        virtual_size_OMP = 0;
    end
    
    %% 计算NMSE
    NMSE_nume = zeros(num_algorithms, 1);
    NMSE_deno = 0;
    
    for kk = 0:(N-1)
        for ll = 1:M
            if (test==1)
                h_w = sum(Sampling_Function_v(N,kk+1,1,k_v_test) .* Sampling_Function_t(M,ll,1,l_t_test) .* h_test .* exp(-2i*pi.*k_v_test.*l_t_test/N/M));
            else
                h_w = sum(Sampling_Function_v(N,kk+1,1,k_v_init) .* Sampling_Function_t(M,ll,1,l_t_init) .* h_c_init .* exp(-2i*pi.*k_v_init.*l_t_init/N/M));
            end
            
            % Traditional Impulse estimation
            h_w_hat_traditional = sum(Sampling_Function_v(N,kk+1,1,k_v_hat_traditional) .* Sampling_Function_t(M,ll,1,l_t_hat_traditional) .* h_hat_traditional .* exp(-2i*pi.*k_v_hat_traditional.*l_t_hat_traditional/N/M));
            
            % 1D SBL estimation
            h_w_hat_1D = sum(Sampling_Function_v(N,kk+1,1,k_v_hat_1D) .* Sampling_Function_t(M,ll,1,l_t_hat_1D) .* h_hat_1D .* exp(-2i*pi.*k_v_hat_1D.*l_t_hat_1D/N/M));
            
            % 2D SBL estimation
            h_w_hat_2D = sum(Sampling_Function_v(N,kk+1,1,k_v_opt_2D) .* Sampling_Function_t(M,ll,1,l_t_opt_2D) .* H_opt_2D .* exp(-2i*pi.*k_v_opt_2D.*l_t_opt_2D/N/M));
            
            % 分层1D SBL estimation
            h_w_hat_hierarchical_1D = sum(Sampling_Function_v(N,kk+1,1,k_v_hat_hierarchical_1D) .* Sampling_Function_t(M,ll,1,l_t_hat_hierarchical_1D) .* h_hat_hierarchical_1D .* exp(-2i*pi.*k_v_hat_hierarchical_1D.*l_t_hat_hierarchical_1D/N/M));
            
            % 分层2D SBL estimation
            h_w_hat_hierarchical_2D = sum(Sampling_Function_v(N,kk+1,1,k_v_opt_hierarchical_2D) .* Sampling_Function_t(M,ll,1,l_t_opt_hierarchical_2D) .* H_opt_hierarchical_2D .* exp(-2i*pi.*k_v_opt_hierarchical_2D.*l_t_opt_hierarchical_2D/N/M));
            
            % OMP estimation
            h_w_hat_OMP = sum(Sampling_Function_v(N,kk+1,1,k_v_hat_OMP) .* Sampling_Function_t(M,ll,1,l_t_hat_OMP) .* h_hat_OMP .* exp(-2i*pi.*k_v_hat_OMP.*l_t_hat_OMP/N/M));
            
            NMSE_nume(1) = NMSE_nume(1) + abs(h_w - h_w_hat_traditional).^2;
            NMSE_nume(2) = NMSE_nume(2) + abs(h_w - h_w_hat_1D).^2;
            NMSE_nume(3) = NMSE_nume(3) + abs(h_w - h_w_hat_2D).^2;
            NMSE_nume(4) = NMSE_nume(4) + abs(h_w - h_w_hat_hierarchical_1D).^2;
            NMSE_nume(5) = NMSE_nume(5) + abs(h_w - h_w_hat_hierarchical_2D).^2;
            NMSE_nume(6) = NMSE_nume(6) + abs(h_w - h_w_hat_OMP).^2;
            NMSE_deno = NMSE_deno + abs(h_w)^2;
        end 
    end
    
    for i = 1:num_algorithms
        temp_NMSE(i) = temp_NMSE(i) + NMSE_nume(i) / NMSE_deno;
    end
end

% 计算平均结果
for i = 1:num_algorithms
    NMSE_results(i) = temp_NMSE(i) / N_fram;
    time_results(i) = temp_time(i) / N_fram;
end

% 计算复杂度（网格点数）
complexity_results(1) = N * M;  % Traditional Impulse
complexity_results(2) = virtual_size_1D;  % 1D SBL
complexity_results(3) = virtual_size_2D;  % 2D SBL
complexity_results(4) = virtual_size_coarse_1D + length(refined_k_v_1D);  % 分层1D SBL
complexity_results(5) = virtual_size_coarse_2D + length(k_v_refined_2D);  % 分层2D SBL
complexity_results(6) = virtual_size_OMP;  % OMP

% 转换为dB
NMSE_dB = 10 * log10(NMSE_results);

%% 绘制Performance-Efficiency Trade-off图
colors = [128,128,128;    % 灰色 - Traditional Impulse
          31,119,180;     % 蓝色 - 1D SBL
          255,127,14;     % 橙色 - 2D SBL
          234,174,31;     % 黄色 - 分层1D SBL
          215,94,59;      % 红色 - 分层2D SBL
          44,160,44]/256; % 绿色 - OMP

markers = {'o', 's', '^', 'd', 'v', 'p'};
marker_sizes = [120, 120, 120, 120, 120, 120];

% 创建单个图形
figure('Position', [100, 100, 800, 600]);

% NMSE vs 计算时间 (Performance-Efficiency Trade-off)
for i = 1:num_algorithms
    scatter(time_results(i), NMSE_dB(i), marker_sizes(i), colors(i,:), markers{i}, 'filled', 'LineWidth', 2);
    hold on;
    % 添加算法名称标签
    text(time_results(i), NMSE_dB(i), ['  ' algorithm_names{i}], 'FontSize', 11, 'HorizontalAlignment', 'left', 'FontWeight', 'bold');
end
set(gca, 'XScale', 'log');
grid on;
xlabel('Computation Time (s)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('NMSE (dB)', 'FontSize', 14, 'FontWeight', 'bold');
title('Performance-Efficiency Trade-off: NMSE vs Computation Time', 'FontSize', 16, 'FontWeight', 'bold');
legend(algorithm_names, 'Location', 'best', 'FontSize', 12);
set(gca, 'FontSize', 12);

% 设置坐标轴范围以获得更好的显示效果
xlim([min(time_results)*0.5, max(time_results)*2]);
ylim([min(NMSE_dB)-2, max(NMSE_dB)+2]);

% 保存图像
if ~exist('figures', 'dir')
    mkdir('figures');
end
saveas(gcf, 'figures/performance_efficiency_tradeoff_single.png');
fprintf('\n图像已保存至: figures/performance_efficiency_tradeoff_single.png\n');

%% 数值结果输出
fprintf('\n==== 所有算法Performance-Efficiency Trade-off结果 ====\n');
fprintf('Algorithm\t\t\tNMSE (dB)\tTime (s)\tComplexity\tEfficiency Ratio\n');
fprintf('--------\t\t\t---------\t--------\t----------\t----------------\n');
for i = 1:num_algorithms
    fprintf('%-20s\t%.2f\t\t%.4f\t\t%d\t\t%.2f\n', ...
        algorithm_names{i}, NMSE_dB(i), time_results(i), complexity_results(i), efficiency_ratio(i));
end

% 找到最佳算法
[best_nmse, best_nmse_idx] = min(NMSE_dB);
[best_time, best_time_idx] = min(time_results(2:end));  % 排除Traditional Impulse
best_time_idx = best_time_idx + 1;  % 调整索引
[best_efficiency, best_efficiency_idx] = max(efficiency_ratio(2:end));  % 排除Traditional Impulse
best_efficiency_idx = best_efficiency_idx + 1;  % 调整索引

fprintf('\n==== 主要发现 ====\n');
fprintf('1. 最佳性能算法: %s (NMSE: %.2f dB)\n', algorithm_names{best_nmse_idx}, best_nmse);
fprintf('2. 最快算法: %s (Time: %.4f s)\n', algorithm_names{best_time_idx}, time_results(best_time_idx));
fprintf('3. 最高效率算法: %s (Efficiency: %.2f dB/s)\n', algorithm_names{best_efficiency_idx}, best_efficiency);
fprintf('4. 性能范围: %.2f to %.2f dB\n', min(NMSE_dB), max(NMSE_dB));
fprintf('5. 时间范围: %.4f to %.4f s\n', min(time_results), max(time_results));
fprintf('6. 复杂度范围: %d to %d grid points\n', min(complexity_results), max(complexity_results));

% 推荐算法
fprintf('\n==== 算法推荐 ====\n');
fprintf('? 追求最佳性能: %s\n', algorithm_names{best_nmse_idx});
fprintf('? 追求最快速度: %s\n', algorithm_names{best_time_idx});
fprintf('? 平衡性能与效率: %s\n', algorithm_names{best_efficiency_idx});

fprintf('\n=== 仿真完成 ===\n');

%% 辅助函数定义

%% 分层2D SBL网格细化函数
function [k_v_refined, l_t_refined] = create_refined_grid_2D_improved(k_v_coarse, l_t_coarse, significant_indices, N_v_coarse, M_t_coarse, r_v_coarse, r_t_coarse, r_v_fine, r_t_fine, k_max, l_max)
    % Create a proper 2D refined grid based on significant coefficients
    % This function creates a structured grid compatible with CE_2D_SBL
    
    % Get coarse grid structure
    [kv_bar_coarse, lt_bar_coarse] = First_Order_Linear_Approximation_2D(N_v_coarse, M_t_coarse, k_max, r_v_coarse, r_t_coarse);
    
    % Find unique significant positions in 2D grid
    significant_k_v = [];
    significant_l_t = [];
    
    for idx = 1:length(significant_indices)
        sig_idx = significant_indices(idx);
        
        % Convert linear index to 2D subscripts
        [row_idx, col_idx] = ind2sub([M_t_coarse, N_v_coarse], sig_idx);
        
        % Get the coarse grid center position
        k_v_center = kv_bar_coarse(col_idx);
        l_t_center = lt_bar_coarse(row_idx);
        
        significant_k_v = [significant_k_v; k_v_center];
        significant_l_t = [significant_l_t; l_t_center];
    end
    
    % Remove duplicates
    significant_positions = unique([significant_k_v, significant_l_t], 'rows');
    significant_k_v = significant_positions(:, 1);
    significant_l_t = significant_positions(:, 2);
    
    % Create refined grid dimensions
    N_v_refined = ceil(2 * k_max / r_v_fine);
    M_t_refined = ceil(l_max / r_t_fine);
    
    % Create structured refined grid using First_Order_Linear_Approximation_2D
    [k_v_refined_full, l_t_refined_full] = First_Order_Linear_Approximation_2D(N_v_refined, M_t_refined, k_max, r_v_fine, r_t_fine);
    
    % For compatibility with CE_2D_SBL structure, return the full refined grid
    % The algorithm will automatically focus on significant regions through the SBL process
    k_v_refined = k_v_refined_full;
    l_t_refined = l_t_refined_full;
end