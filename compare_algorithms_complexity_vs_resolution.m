close all
clear all
clc

%% 算法复杂度随网格分辨率变化分析
% 固定SNR=10dB，测试不同网格分辨率下各算法的计算时间
% 包含：1D SBL, 2D SBL, 分层1D SBL, 分层2D SBL, OMP, Traditional Impulse

%% 网格分辨率参数设置
% 测试不同的分辨率值
resolution_values = [0.6, 0.4, 0.3, 0.2, 0.1];  % 从粗到细的分辨率
num_resolutions = length(resolution_values);

% 分层算法的粗细网格比例
coarse_to_fine_ratio = 2.5;  % 粗网格分辨率 = 细网格分辨率 * ratio
threshold_ratio = 0.1;       % 分层算法系数选择阈值

%% 固定参数设置
% 固定SNR
SNR_dB_fixed = 10;  % 固定SNR为10dB
SNR_fixed = 10^(SNR_dB_fixed/10);
sigma_2_fixed = 0.5 / SNR_fixed;
sigma_2_p_fixed = sigma_2_fixed / 1e5;

%% Test Mode
test = 0;    % set to 1 when testing

%% OTFS parameters
N = 32;      % number of symbols in time
M = 32;      % number of subcarriers in frequency
M_mod = 4;   % size of QAM constellation
M_bits = log2(M_mod);
eng_sqrt = (M_mod==2)+(M_mod~=2)*sqrt((M_mod-1)/6*(2^2));

%% delay-Doppler grid symbol placement
car_fre = 5*10^9;  % Carrier frequency
delta_f = 15*10^3; % subcarrier spacing: 15 KHz
T = 1/delta_f;     % one time symbol duration in OTFS frame
k_max = 4;         % system setting
l_max = 4;
t_max = l_max / M / delta_f;
v_max = k_max / N / T;
P = 9;
on_flag = 0;      % set to 0 when on-grid

%% pilot symbol placing
Kp = 1;
Lp = 1;
x_kp = floor(N/2);
x_lp = floor(M/2);
data_grid = ones(N,M);
data_grid(x_kp-floor(Kp/2)-2*k_max:x_kp+floor(Kp/2)+2*k_max,x_lp-floor(Lp/2)-l_max:x_lp+floor(Lp/2)+l_max)=0;
N_syms_perfram = sum(sum(data_grid));
N_bits_perfram = N_syms_perfram*M_bits;

%% 仿真参数
N_fram = 5;  % 减少仿真帧数以专注于复杂度测试

%% 初始化结果存储
% 计算时间存储
time_1D_SBL = zeros(num_resolutions, 1);
time_2D_SBL = zeros(num_resolutions, 1);
time_hierarchical_1D = zeros(num_resolutions, 1);
time_hierarchical_2D = zeros(num_resolutions, 1);
time_OMP = zeros(num_resolutions, 1);
time_traditional = zeros(num_resolutions, 1);

% 网格点数存储
grid_points_1D = zeros(num_resolutions, 1);
grid_points_2D = zeros(num_resolutions, 1);
grid_points_hierarchical_1D_coarse = zeros(num_resolutions, 1);
grid_points_hierarchical_1D_fine = zeros(num_resolutions, 1);
grid_points_hierarchical_2D_coarse = zeros(num_resolutions, 1);
grid_points_hierarchical_2D_fine = zeros(num_resolutions, 1);
grid_points_OMP = zeros(num_resolutions, 1);
grid_points_traditional = zeros(num_resolutions, 1);

fprintf('\n=== 算法复杂度随网格分辨率变化分析 ===\n');
fprintf('固定SNR: %d dB\n', SNR_dB_fixed);
fprintf('仿真参数：N=%d, M=%d, P=%d, N_fram=%d\n', N, M, P, N_fram);
fprintf('测试分辨率范围: %.2f ~ %.2f\n', max(resolution_values), min(resolution_values));
fprintf('\n');

%% 主循环：遍历不同分辨率
for res_idx = 1:num_resolutions
    r_current = resolution_values(res_idx);
    
    % 设置当前分辨率下的参数
    r_v_1D = r_current;
    r_t_1D = r_current;
    r_v_2D = r_current;
    r_t_2D = r_current;
    r_v_OMP = r_current;
    r_t_OMP = r_current;
    
    % 分层算法参数
    r_v_coarse_1D = r_current * coarse_to_fine_ratio;
    r_t_coarse_1D = r_current * coarse_to_fine_ratio;
    r_v_fine_1D = r_current;
    r_t_fine_1D = r_current;
    
    r_v_coarse_2D = r_current * coarse_to_fine_ratio;
    r_t_coarse_2D = r_current * coarse_to_fine_ratio;
    r_v_fine_2D = r_current;
    r_t_fine_2D = r_current;
    
    fprintf('\n--- 测试分辨率: %.3f ---\n', r_current);
    
    % 计算网格点数
    grid_points_1D(res_idx) = ceil(2 * k_max / r_v_1D) * ceil(l_max / r_t_1D);
    grid_points_2D(res_idx) = ceil(2 * k_max / r_v_2D) * ceil(l_max / r_t_2D);
    grid_points_hierarchical_1D_coarse(res_idx) = ceil(2 * k_max / r_v_coarse_1D) * ceil(l_max / r_t_coarse_1D);
    grid_points_hierarchical_1D_fine(res_idx) = ceil(2 * k_max / r_v_fine_1D) * ceil(l_max / r_t_fine_1D);
    grid_points_hierarchical_2D_coarse(res_idx) = ceil(2 * k_max / r_v_coarse_2D) * ceil(l_max / r_t_coarse_2D);
    grid_points_hierarchical_2D_fine(res_idx) = ceil(2 * k_max / r_v_fine_2D) * ceil(l_max / r_t_fine_2D);
    grid_points_OMP(res_idx) = ceil(2 * k_max / r_v_OMP) * ceil(l_max / r_t_OMP);
    grid_points_traditional(res_idx) = (2 * k_max + 1) * (l_max + 1);
    
    fprintf('网格点数 - 1D SBL: %d, 2D SBL: %d, OMP: %d\n', ...
        grid_points_1D(res_idx), grid_points_2D(res_idx), grid_points_OMP(res_idx));
    
    % 初始化时间累计器
    total_time_1D = 0;
    total_time_2D = 0;
    total_time_hierarchical_1D = 0;
    total_time_hierarchical_2D = 0;
    total_time_OMP = 0;
    total_time_traditional = 0;
    
    %% Monte Carlo仿真
    for ifram = 1:N_fram
        if mod(ifram, 2) == 0
            fprintf('  Frame %d/%d\n', ifram, N_fram);
        end
        
        %% 随机信道初始化
        v_c_init = unifrnd(-v_max,v_max,P,1);
        t_c_init = unifrnd(0,t_max,P,1);
        l_ti = t_c_init.*(M * delta_f);
        q_l_t = exp(-0.1.*l_ti)./sum(exp(-0.1.*l_ti));
        h_c_init = normrnd(0,q_l_t);
        k_v_init = v_c_init .*(N*T);
        l_t_init = t_c_init .*(M*delta_f);
        
        pow_prof = (1/P) * (ones(1,P));
        chan_coef = zeros(1,P);
        
        %% 随机输入比特生成
        trans_info_bit = randi([0,1],N_syms_perfram*M_bits,1);
        data=qammod(reshape(trans_info_bit,M_bits,N_syms_perfram), M_mod,'gray','InputType','bit');
        data = data./(eng_sqrt);
        x = Generate_2D_data_grid_CE(N,M,data,data_grid);
        x(x_kp-floor(Kp/2):x_kp+floor(Kp/2),x_lp-floor(Lp/2):x_lp+floor(Lp/2))=sqrt(1000)/(Kp*Lp);
        
        % 计算测量矩阵
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
        
        %% 信道输出
        noise_gen_re = sigma_2_fixed * randn(M*N,1);
        noise_gen_im = sigma_2_fixed * randn(M*N,1);
        noise_gen = noise_gen_re + 1i.* noise_gen_im;
        h_c_init_hat = h_c_init .* exp(-2i * pi/M/N * (k_v_init.*l_t_init));
        r = phi_sys * h_c_init + noise_gen;
        y = reshape(r,M,N).';
        
        %% 获取截断测量矩阵
        y_trunc = y(floor(N/2)-floor(Kp/2)-k_max:floor(N/2)+floor(Kp/2)+k_max,floor(M/2)-floor(Lp/2):floor(M/2)+floor(Lp/2)+l_max);
        N_T = 2*k_max+Kp;
        M_T = l_max+Lp;
        y_T = reshape(y_trunc.',N_T*M_T,1);
        
        %% 算法1: 1D SBL Channel Estimation
        tic;
        try
            [h_hat_1D,k_v_hat_1D,l_t_hat_1D,virtual_size_1D, Phi_1D, delta_a] = CE_1D_SBL(sqrt(1000)/(Kp*Lp),k_max+Kp,Lp,M,N,N_T,M_T,y_T,r_v_1D,r_t_1D,k_max,l_max,on_flag);
        catch ME
            fprintf('    1D SBL算法出错: %s\n', ME.message);
        end
        total_time_1D = total_time_1D + toc;
        
        %% 算法2: 2D SBL Channel Estimation
        tic;
        try
            [H_opt_2D,k_v_opt_2D,l_t_opt_2D,virtual_size_2D] = CE_2D_SBL(sqrt(1000)/(2*Kp*Lp),k_max+Kp,2,Lp,M,N,N_T,M_T,y_trunc,r_v_2D,r_t_2D,k_max,l_max,on_flag);
        catch ME
            fprintf('    2D SBL算法出错: %s\n', ME.message);
        end
        total_time_2D = total_time_2D + toc;
        
        %% 算法3: 分层1D SBL Channel Estimation
        tic;
        try
            % Step 1: Coarse Grid SBL Estimation
            [h_hat_coarse_1D, k_v_hat_coarse_1D, l_t_hat_coarse_1D, virtual_size_coarse_1D, Phi_coarse_1D, delta_coarse_1D] = ...
                CE_1D_SBL(sqrt(1000)/(Kp*Lp), k_max+Kp, Lp, M, N, N_T, M_T, y_T, r_v_coarse_1D, r_t_coarse_1D, k_max, l_max, on_flag);
            
            % Step 2: Identify significant coefficients and their locations
            [sorted_coeff_1D, sorted_idx_1D] = sort(abs(h_hat_coarse_1D), 'descend');
            num_significant_1D = max(1, floor(threshold_ratio * length(h_hat_coarse_1D)));
            significant_idx_1D = sorted_idx_1D(1:num_significant_1D);
            
            % Extract corresponding k_v and l_t values for significant coefficients
            k_v_significant_1D = k_v_hat_coarse_1D(significant_idx_1D);
            l_t_significant_1D = l_t_hat_coarse_1D(significant_idx_1D);
            
            % Step 3: Create refined grid around significant locations
            refined_k_v_1D = [];
            refined_l_t_1D = [];
            
            for i = 1:length(k_v_significant_1D)
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
        catch ME
            fprintf('    分层1D SBL算法出错: %s\n', ME.message);
        end
        total_time_hierarchical_1D = total_time_hierarchical_1D + toc;
        
        %% 算法4: 分层2D SBL Channel Estimation
        tic;
        try
            % Step 1: Coarse Grid 2D SBL Estimation
            [H_opt_coarse_2D,k_v_opt_coarse_2D,l_t_opt_coarse_2D,virtual_size_coarse_2D] = CE_2D_SBL(sqrt(1000)/(2*Kp*Lp),k_max+Kp,2,Lp,M,N,N_T,M_T,y_trunc,r_v_coarse_2D,r_t_coarse_2D,k_max,l_max,on_flag);
            
            % Step 2: Identify significant coefficients and their 2D positions
            coeff_magnitude_2D = abs(H_opt_coarse_2D);
            max_coeff_2D = max(coeff_magnitude_2D);
            significant_indices_2D = find(coeff_magnitude_2D > threshold_ratio * max_coeff_2D);
            
            % Convert linear indices to 2D grid positions
            N_v_coarse_2D = ceil(2 * k_max / r_v_coarse_2D);
            M_t_coarse_2D = ceil(l_max / r_t_coarse_2D);
            [row_indices_2D, col_indices_2D] = ind2sub([M_t_coarse_2D, N_v_coarse_2D], significant_indices_2D);
            
            % Step 3: Create refined grid around significant coefficients
            [k_v_refined_2D, l_t_refined_2D] = create_refined_grid_2D_improved(k_v_opt_coarse_2D, l_t_opt_coarse_2D, significant_indices_2D, N_v_coarse_2D, M_t_coarse_2D, r_v_coarse_2D, r_t_coarse_2D, r_v_fine_2D, r_t_fine_2D, k_max, l_max);
            
            % Step 4: Refined 2D SBL estimation on the refined grid
            [H_opt_hierarchical_2D, k_v_opt_hierarchical_2D, l_t_opt_hierarchical_2D] = hierarchical_SBL_refined_2D(sqrt(1000)/(2*Kp*Lp), k_max+Kp, 2, Lp, M, N, N_T, M_T, y_trunc, k_v_refined_2D, l_t_refined_2D, r_v_fine_2D, r_t_fine_2D, k_max, l_max, on_flag);
        catch ME
            fprintf('    分层2D SBL算法出错: %s\n', ME.message);
        end
        total_time_hierarchical_2D = total_time_hierarchical_2D + toc;
        
        %% 算法5: OMP Channel Estimation
        tic;
        try
            [h_hat_OMP, omp_index, k_v_hat_OMP, l_t_hat_OMP] = OMP(sqrt(1000)/(Kp*Lp),k_max+Kp,Lp,M,N,N_T,M_T,y_T,r_v_OMP,r_t_OMP,k_max,l_max);
        catch ME
            fprintf('    OMP算法出错: %s\n', ME.message);
        end
        total_time_OMP = total_time_OMP + toc;
        
        %% 算法6: Traditional Impulse Channel Estimation
        tic;
        try
            [h_hat_traditional, k_v_hat_traditional, l_t_hat_traditional] = traditional_impulse(sqrt(1000)/(Kp*Lp), y_trunc, k_max, l_max, sigma_2_p_fixed);
        catch ME
            fprintf('    Traditional算法出错: %s\n', ME.message);
        end
        total_time_traditional = total_time_traditional + toc;
    end
    
    % 计算平均时间
    time_1D_SBL(res_idx) = total_time_1D / N_fram;
    time_2D_SBL(res_idx) = total_time_2D / N_fram;
    time_hierarchical_1D(res_idx) = total_time_hierarchical_1D / N_fram;
    time_hierarchical_2D(res_idx) = total_time_hierarchical_2D / N_fram;
    time_OMP(res_idx) = total_time_OMP / N_fram;
    time_traditional(res_idx) = total_time_traditional / N_fram;
    
    fprintf('平均计算时间 (秒): 1D SBL=%.4f, 2D SBL=%.4f, 分层1D=%.4f, 分层2D=%.4f, OMP=%.4f, Traditional=%.4f\n', ...
        time_1D_SBL(res_idx), time_2D_SBL(res_idx), time_hierarchical_1D(res_idx), ...
        time_hierarchical_2D(res_idx), time_OMP(res_idx), time_traditional(res_idx));
end

%% 结果绘制
colors = [0,107,182; %蓝1
          118,41,133; %紫2
          234,174,31; %黄3
          215,94,59; %橙4
          184,125,182; %粉紫5
          71,90,40; %绿6
          161,27,30]; %红7
colors = colors/256;

% 图1: 计算时间随分辨率变化
figure('Position', [100, 100, 1200, 800]);

% 子图1: 计算时间对比
subplot(2, 2, 1);
semilogy(resolution_values, time_1D_SBL, '-v', 'Color', colors(3,:), 'LineWidth', 2, 'MarkerSize', 8);
hold on;
semilogy(resolution_values, time_2D_SBL, '-*', 'Color', colors(2,:), 'LineWidth', 2, 'MarkerSize', 8);
semilogy(resolution_values, time_hierarchical_1D, '-s', 'Color', colors(1,:), 'LineWidth', 2, 'MarkerSize', 8);
semilogy(resolution_values, time_hierarchical_2D, '-o', 'Color', colors(4,:), 'LineWidth', 2, 'MarkerSize', 8);
semilogy(resolution_values, time_OMP, '-d', 'Color', colors(6,:), 'LineWidth', 2, 'MarkerSize', 8);
semilogy(resolution_values, time_traditional, '-x', 'Color', colors(7,:), 'LineWidth', 2, 'MarkerSize', 8);
grid on;
xlabel('网格分辨率');
ylabel('平均计算时间 (秒)');
title('算法计算时间随网格分辨率变化');
legend('1D SBL', '2D SBL', '分层1D SBL', '分层2D SBL', 'OMP', 'Traditional', 'Location', 'best');
set(gca, 'XDir', 'reverse');  % 反转x轴，使分辨率从粗到细

% 子图2: 网格点数对比
subplot(2, 2, 2);
semilogy(resolution_values, grid_points_1D, '-v', 'Color', colors(3,:), 'LineWidth', 2, 'MarkerSize', 8);
hold on;
semilogy(resolution_values, grid_points_2D, '-*', 'Color', colors(2,:), 'LineWidth', 2, 'MarkerSize', 8);
semilogy(resolution_values, grid_points_OMP, '-d', 'Color', colors(6,:), 'LineWidth', 2, 'MarkerSize', 8);
semilogy(resolution_values, grid_points_traditional, '-x', 'Color', colors(7,:), 'LineWidth', 2, 'MarkerSize', 8);
grid on;
xlabel('网格分辨率');
ylabel('网格点数');
title('网格点数随分辨率变化');
legend('1D SBL', '2D SBL', 'OMP', 'Traditional', 'Location', 'best');
set(gca, 'XDir', 'reverse');

% 子图3: 计算效率对比 (时间/网格点数)
subplot(2, 2, 3);
efficiency_1D = time_1D_SBL ./ grid_points_1D * 1e6;  % 微秒/网格点
efficiency_2D = time_2D_SBL ./ grid_points_2D * 1e6;
efficiency_OMP = time_OMP ./ grid_points_OMP * 1e6;
efficiency_traditional = time_traditional ./ grid_points_traditional * 1e6;

plot(resolution_values, efficiency_1D, '-v', 'Color', colors(3,:), 'LineWidth', 2, 'MarkerSize', 8);
hold on;
plot(resolution_values, efficiency_2D, '-*', 'Color', colors(2,:), 'LineWidth', 2, 'MarkerSize', 8);
plot(resolution_values, efficiency_OMP, '-d', 'Color', colors(6,:), 'LineWidth', 2, 'MarkerSize', 8);
plot(resolution_values, efficiency_traditional, '-x', 'Color', colors(7,:), 'LineWidth', 2, 'MarkerSize', 8);
grid on;
xlabel('网格分辨率');
ylabel('计算效率 (微秒/网格点)');
title('算法计算效率对比');
legend('1D SBL', '2D SBL', 'OMP', 'Traditional', 'Location', 'best');
set(gca, 'XDir', 'reverse');

% 子图4: 分层算法粗细网格点数对比
subplot(2, 2, 4);
semilogy(resolution_values, grid_points_hierarchical_1D_coarse, '--s', 'Color', colors(1,:), 'LineWidth', 2, 'MarkerSize', 6);
hold on;
semilogy(resolution_values, grid_points_hierarchical_1D_fine, '-s', 'Color', colors(1,:), 'LineWidth', 2, 'MarkerSize', 8);
semilogy(resolution_values, grid_points_hierarchical_2D_coarse, '--o', 'Color', colors(4,:), 'LineWidth', 2, 'MarkerSize', 6);
semilogy(resolution_values, grid_points_hierarchical_2D_fine, '-o', 'Color', colors(4,:), 'LineWidth', 2, 'MarkerSize', 8);
grid on;
xlabel('网格分辨率');
ylabel('网格点数');
title('分层算法粗细网格点数');
legend('分层1D粗网格', '分层1D细网格', '分层2D粗网格', '分层2D细网格', 'Location', 'best');
set(gca, 'XDir', 'reverse');

sgtitle(sprintf('算法复杂度随网格分辨率变化分析 (SNR=%ddB)', SNR_dB_fixed), 'FontSize', 14, 'FontWeight', 'bold');

%% 数值结果表格
fprintf('\n=== 算法复杂度随网格分辨率变化结果 ===\n');
fprintf('分辨率\t网格点数\t\t\t\t\t\t\t\t计算时间(秒)\n');
fprintf('\t\t1D SBL\t2D SBL\tOMP\t\tTraditional\t1D SBL\t\t2D SBL\t\t分层1D\t\t分层2D\t\tOMP\t\t\tTraditional\n');
fprintf('------\t------\t------\t------\t----------\t------\t\t------\t\t------\t\t------\t\t------\t\t----------\n');
for i = 1:num_resolutions
    fprintf('%.3f\t%d\t\t%d\t\t%d\t\t%d\t\t\t%.4f\t\t%.4f\t\t%.4f\t\t%.4f\t\t%.4f\t\t%.4f\n', ...
        resolution_values(i), grid_points_1D(i), grid_points_2D(i), grid_points_OMP(i), grid_points_traditional(i), ...
        time_1D_SBL(i), time_2D_SBL(i), time_hierarchical_1D(i), time_hierarchical_2D(i), time_OMP(i), time_traditional(i));
end

%% 复杂度增长率分析
fprintf('\n=== 复杂度增长率分析 ===\n');
for i = 2:num_resolutions
    ratio_resolution = resolution_values(i-1) / resolution_values(i);
    ratio_time_1D = time_1D_SBL(i) / time_1D_SBL(i-1);
    ratio_time_2D = time_2D_SBL(i) / time_2D_SBL(i-1);
    ratio_grid_1D = grid_points_1D(i) / grid_points_1D(i-1);
    ratio_grid_2D = grid_points_2D(i) / grid_points_2D(i-1);
    
    fprintf('分辨率提升%.2fx: 1D SBL时间增长%.2fx(网格点%.2fx), 2D SBL时间增长%.2fx(网格点%.2fx)\n', ...
        ratio_resolution, ratio_time_1D, ratio_grid_1D, ratio_time_2D, ratio_grid_2D);
end

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