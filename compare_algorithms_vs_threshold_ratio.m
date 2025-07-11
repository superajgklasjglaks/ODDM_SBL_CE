close all;
clear all;
clc;

%% 分层SBL算法threshold_ratio参数分析
% 基于compare_all_algorithms.m代码结构，只保留分层1D和分层2D SBL算法

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

%% threshold_ratio参数范围
threshold_ratio_values = 0.05:0.05:0.25;  % 从0.05到0.25，步长0.05
num_threshold = length(threshold_ratio_values);

%% 分层算法分辨率参数设置
% 分层1D SBL 参数
r_v_coarse_1D = 0.5;    % 分层1D SBL 粗网格多普勒分辨率
r_t_coarse_1D = 0.5;    % 分层1D SBL 粗网格时延分辨率
r_v_fine_1D = 0.2;      % 分层1D SBL 细网格多普勒分辨率
r_t_fine_1D = 0.2;      % 分层1D SBL 细网格时延分辨率

% 分层2D SBL 参数
r_v_coarse_2D = 0.5;    % 分层2D SBL 粗网格多普勒分辨率
r_t_coarse_2D = 0.5;    % 分层2D SBL 粗网格时延分辨率
r_v_fine_2D = 0.2;      % 分层2D SBL 细网格多普勒分辨率
r_t_fine_2D = 0.2;      % 分层2D SBL 细网格时延分辨率

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
N_fram = 10;  % 仿真帧数
NMSE_count_hierarchical_1D = zeros(num_threshold,1);
NMSE_count_hierarchical_2D = zeros(num_threshold,1);
time_hierarchical_1D = zeros(num_threshold,1);
time_hierarchical_2D = zeros(num_threshold,1);

fprintf('\n=== 分层SBL算法threshold_ratio参数分析开始 ===\n');
fprintf('仿真参数：N=%d, M=%d, P=%d, N_fram=%d, SNR=%d dB\n', N, M, P, N_fram, SNR_dB_fixed);
fprintf('threshold_ratio范围: %.2f - %.2f (步长: %.2f)\n', min(threshold_ratio_values), max(threshold_ratio_values), threshold_ratio_values(2)-threshold_ratio_values(1));
fprintf('\n');

for ithreshold = 1:num_threshold
    threshold_ratio_1D = threshold_ratio_values(ithreshold);
    threshold_ratio_2D = threshold_ratio_values(ithreshold);
    
    fprintf('\n正在仿真 threshold_ratio = %.2f (%d/%d)\n', threshold_ratio_1D, ithreshold, num_threshold);
    
    % 临时存储每帧仿真结果
    temp_NMSE_1D = 0;
    temp_NMSE_2D = 0;
    temp_time_1D = 0;
    temp_time_2D = 0;
    
    for ifram = 1:N_fram
        if mod(ifram, 5) == 0
            fprintf('  Frame %d/%d\n', ifram, N_fram);
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
        
        %% 算法1: 分层1D SBL Channel Estimation
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
            
            temp_time_1D = temp_time_1D + toc;
            
        catch ME
            fprintf('    分层1D SBL算法出错: %s\n', ME.message);
            temp_time_1D = temp_time_1D + toc;
            % 使用默认值
            h_hat_hierarchical_1D = zeros(P,1);
            k_v_hat_hierarchical_1D = zeros(P,1);
            l_t_hat_hierarchical_1D = zeros(P,1);
        end

        %% 算法2: 分层2D SBL Channel Estimation
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
            
            temp_time_2D = temp_time_2D + toc;
            
        catch ME
            fprintf('    分层2D SBL算法出错: %s\n', ME.message);
            temp_time_2D = temp_time_2D + toc;
            % 使用默认值
            H_opt_hierarchical_2D = zeros(P,1);
            k_v_opt_hierarchical_2D = zeros(P,1);
            l_t_opt_hierarchical_2D = zeros(P,1);
        end
        
        %% 计算NMSE
        NMSE_nume_hierarchical_1D = 0;
        NMSE_nume_hierarchical_2D = 0;
        NMSE_deno = 0;
        
        for kk = 0:(N-1)
            for ll = 1:M
                if (test==1)
                    h_w = sum(Sampling_Function_v(N,kk+1,1,k_v_test) .* Sampling_Function_t(M,ll,1,l_t_test) .* h_test .* exp(-2i*pi.*k_v_test.*l_t_test/N/M));
                else
                    h_w = sum(Sampling_Function_v(N,kk+1,1,k_v_init) .* Sampling_Function_t(M,ll,1,l_t_init) .* h_c_init .* exp(-2i*pi.*k_v_init.*l_t_init/N/M));
                end
                
                % 分层1D SBL estimation
                h_w_hat_hierarchical_1D = sum(Sampling_Function_v(N,kk+1,1,k_v_hat_hierarchical_1D) .* Sampling_Function_t(M,ll,1,l_t_hat_hierarchical_1D) .* h_hat_hierarchical_1D .* exp(-2i*pi.*k_v_hat_hierarchical_1D.*l_t_hat_hierarchical_1D/N/M));
                
                % 分层2D SBL estimation
                h_w_hat_hierarchical_2D = sum(Sampling_Function_v(N,kk+1,1,k_v_opt_hierarchical_2D) .* Sampling_Function_t(M,ll,1,l_t_opt_hierarchical_2D) .* H_opt_hierarchical_2D .* exp(-2i*pi.*k_v_opt_hierarchical_2D.*l_t_opt_hierarchical_2D/N/M));
                
                NMSE_nume_hierarchical_1D = NMSE_nume_hierarchical_1D + abs(h_w - h_w_hat_hierarchical_1D).^2;
                NMSE_nume_hierarchical_2D = NMSE_nume_hierarchical_2D + abs(h_w - h_w_hat_hierarchical_2D).^2;
                NMSE_deno = NMSE_deno + abs(h_w)^2;
            end 
        end
        
        temp_NMSE_1D = temp_NMSE_1D + NMSE_nume_hierarchical_1D / NMSE_deno;
        temp_NMSE_2D = temp_NMSE_2D + NMSE_nume_hierarchical_2D / NMSE_deno;
    end
    
    % 计算平均结果
    NMSE_count_hierarchical_1D(ithreshold) = temp_NMSE_1D / N_fram;
    NMSE_count_hierarchical_2D(ithreshold) = temp_NMSE_2D / N_fram;
    time_hierarchical_1D(ithreshold) = temp_time_1D / N_fram;
    time_hierarchical_2D(ithreshold) = temp_time_2D / N_fram;
    
    fprintf('  结果: 分层1D NMSE=%.2f dB, 分层2D NMSE=%.2f dB\n', ...
        10*log10(NMSE_count_hierarchical_1D(ithreshold)), 10*log10(NMSE_count_hierarchical_2D(ithreshold)));
end

%% Results plotting
NMSE_dB_hierarchical_1D = 10 * log10(NMSE_count_hierarchical_1D);
NMSE_dB_hierarchical_2D = 10 * log10(NMSE_count_hierarchical_2D);

colors = [234,174,31;     % 黄色 - 分层1D SBL
          215,94,59]/256; % 橙色 - 分层2D SBL

% 创建包含4个子图的图形
figure('Position', [100, 100, 1200, 800]);

% 子图1: NMSE性能对比
subplot(2, 2, 1);
plot(threshold_ratio_values, NMSE_dB_hierarchical_1D, '-^', 'Color', colors(1,:), 'LineWidth', 2, 'MarkerSize', 8);
hold on;
plot(threshold_ratio_values, NMSE_dB_hierarchical_2D, '-v', 'Color', colors(2,:), 'LineWidth', 2, 'MarkerSize', 8);
grid on;
xlabel('Threshold Ratio', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('NMSE (dB)', 'FontSize', 12, 'FontWeight', 'bold');
title('NMSE Performance vs Threshold Ratio', 'FontSize', 14, 'FontWeight', 'bold');
legend('Hierarchical 1D SBL', 'Hierarchical 2D SBL', 'Location', 'best');
set(gca, 'FontSize', 11);

% 子图2: 计算时间对比
subplot(2, 2, 2);
semilogy(threshold_ratio_values, time_hierarchical_1D, '-^', 'Color', colors(1,:), 'LineWidth', 2, 'MarkerSize', 8);
hold on;
semilogy(threshold_ratio_values, time_hierarchical_2D, '-v', 'Color', colors(2,:), 'LineWidth', 2, 'MarkerSize', 8);
grid on;
xlabel('Threshold Ratio', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Computation Time (s)', 'FontSize', 12, 'FontWeight', 'bold');
title('Computation Time vs Threshold Ratio', 'FontSize', 14, 'FontWeight', 'bold');
legend('Hierarchical 1D SBL', 'Hierarchical 2D SBL', 'Location', 'best');
set(gca, 'FontSize', 11);

% 子图3: 性能效率对比（NMSE改善 vs 时间成本）
subplot(2, 2, 3);
% 计算相对于最大threshold_ratio的性能改善
nmse_improvement_1d = NMSE_dB_hierarchical_1D(end) - NMSE_dB_hierarchical_1D;
nmse_improvement_2d = NMSE_dB_hierarchical_2D(end) - NMSE_dB_hierarchical_2D;

scatter(time_hierarchical_1D, nmse_improvement_1d, 100, colors(1,:), '^', 'filled');
hold on;
scatter(time_hierarchical_2D, nmse_improvement_2d, 100, colors(2,:), 'v', 'filled');

% 添加threshold_ratio标签
for i = 1:length(threshold_ratio_values)
    text(time_hierarchical_1D(i), nmse_improvement_1d(i), sprintf('%.2f', threshold_ratio_values(i)), ...
         'FontSize', 8, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    text(time_hierarchical_2D(i), nmse_improvement_2d(i), sprintf('%.2f', threshold_ratio_values(i)), ...
         'FontSize', 8, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
end

set(gca, 'XScale', 'log');
grid on;
xlabel('Computation Time (s)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('NMSE Improvement (dB)', 'FontSize', 12, 'FontWeight', 'bold');
title('Performance-Efficiency Trade-off', 'FontSize', 14, 'FontWeight', 'bold');
legend('Hierarchical 1D SBL', 'Hierarchical 2D SBL', 'Location', 'best');
set(gca, 'FontSize', 11);

% 子图4: 阈值敏感性分析
subplot(2, 2, 4);
% 计算性能变化率
performance_variation_1d = abs(diff(NMSE_dB_hierarchical_1D));
performance_variation_2d = abs(diff(NMSE_dB_hierarchical_2D));
threshold_mid = threshold_ratio_values(1:end-1) + diff(threshold_ratio_values)/2;

bar(threshold_mid, [performance_variation_1d, performance_variation_2d], 'grouped');
colormap([colors(1,:); colors(2,:)]);
grid on;
xlabel('Threshold Ratio (midpoint)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('NMSE Variation (dB)', 'FontSize', 12, 'FontWeight', 'bold');
title('Threshold Sensitivity Analysis', 'FontSize', 14, 'FontWeight', 'bold');
legend('Hierarchical 1D SBL', 'Hierarchical 2D SBL', 'Location', 'best');
set(gca, 'FontSize', 11);

% 调整子图间距
sgtitle('Hierarchical SBL Algorithms Performance vs Threshold Ratio Analysis', 'FontSize', 16, 'FontWeight', 'bold');

% 保存图像
if ~exist('figures', 'dir')
    mkdir('figures');
end
saveas(gcf, 'figures/threshold_ratio_analysis.png');
fprintf('\n图像已保存至: figures/threshold_ratio_analysis.png\n');

%% 数值结果输出
fprintf('\n==== 数值结果分析 ====\n');
fprintf('Threshold Ratio\t1D SBL NMSE\t2D SBL NMSE\t1D Time\t\t2D Time\n');
fprintf('---------------\t-----------\t-----------\t-------\t\t-------\n');
for i = 1:length(threshold_ratio_values)
    fprintf('%.2f\t\t%.2f dB\t\t%.2f dB\t\t%.4f s\t\t%.4f s\n', ...
        threshold_ratio_values(i), NMSE_dB_hierarchical_1D(i), NMSE_dB_hierarchical_2D(i), ...
        time_hierarchical_1D(i), time_hierarchical_2D(i));
end

% 找到最佳threshold_ratio
[~, best_idx_1d] = min(NMSE_dB_hierarchical_1D);
[~, best_idx_2d] = min(NMSE_dB_hierarchical_2D);

fprintf('\n==== 主要发现 ====\n');
fprintf('1. 分层1D SBL最佳threshold_ratio: %.2f (NMSE: %.2f dB)\n', ...
    threshold_ratio_values(best_idx_1d), NMSE_dB_hierarchical_1D(best_idx_1d));
fprintf('2. 分层2D SBL最佳threshold_ratio: %.2f (NMSE: %.2f dB)\n', ...
    threshold_ratio_values(best_idx_2d), NMSE_dB_hierarchical_2D(best_idx_2d));
fprintf('3. 性能范围 - 1D SBL: %.2f to %.2f dB\n', ...
    min(NMSE_dB_hierarchical_1D), max(NMSE_dB_hierarchical_1D));
fprintf('4. 性能范围 - 2D SBL: %.2f to %.2f dB\n', ...
    min(NMSE_dB_hierarchical_2D), max(NMSE_dB_hierarchical_2D));
fprintf('5. 时间范围 - 1D SBL: %.4f to %.4f s\n', ...
    min(time_hierarchical_1D), max(time_hierarchical_1D));
fprintf('6. 时间范围 - 2D SBL: %.4f to %.4f s\n', ...
    min(time_hierarchical_2D), max(time_hierarchical_2D));

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