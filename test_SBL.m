clc;
clear;
close all;

%% 系统参数设置
array_num = 10;              % 阵元数
snapshot_num = 100;          % 快拍数
source_aoa = [-30, -5, 10, 45];  % 信源真实角度
c = 340;                     % 波速
f = 1000;                    % 频率
lambda = c / f;              % 波长
d = 0.5 * lambda;            % 阵元间距
source_num = length(source_aoa); % 信源数
SNR_range = -10:5:20;        % 信噪比范围（dB）
MC_trials = 50;              % 蒙特卡洛仿真次数

%% 算法参数设置
% 普通SBL参数
params_SBL.reso = 5;         % 固定网格分辨率1°
params_SBL.maxiter = 1000;
params_SBL.tolerance = 1e-4;

% 多分辨率SBL参数
params_MultiRes.coarse_reso = 5;  % 粗网格分辨率5°
params_MultiRes.fine_reso = 1;    % 精网格分辨率1°
params_MultiRes.window_size = 10; % 精网格窗口大小20°
params_MultiRes.maxiter = 1000;
params_MultiRes.tolerance = 1e-4;

%% 存储结果
NMSE_SBL = zeros(length(SNR_range), 1);
NMSE_MultiRes = zeros(length(SNR_range), 1);

%% 主循环：不同SNR
for snr_idx = 1:length(SNR_range)
    SNR = SNR_range(snr_idx);
    fprintf('Processing SNR = %d dB...\n', SNR);
    
    % 临时存储误差
    errors_SBL = zeros(MC_trials, 1);
    errors_MultiRes = zeros(MC_trials, 1);
    
    %% 蒙特卡洛仿真
    for mc = 1:MC_trials
        %% 生成接收信号
        A = exp(-1i * (0:array_num-1)' * 2 * pi * (d / lambda) * sind(source_aoa));
        X = (randn(source_num, snapshot_num) + 1i * randn(source_num, snapshot_num)) / sqrt(2);
        signal_power = mean(abs(A * X).^2, 'all');
        noise_power = signal_power / (10^(SNR / 10));
        n = sqrt(noise_power) * (randn(array_num, snapshot_num) + 1i * randn(array_num, snapshot_num)) / sqrt(2);
        Y = A * X + n;
        
        %% 普通SBL估计
        params_SBL.Y = Y;
        params_SBL.sigma2 = noise_power;
        res_SBL = StandardSBL(params_SBL);
        
        % 提取角度估计
        est_aoa_SBL = extract_angles(res_SBL, source_num);
        errors_SBL(mc) = calculate_error(est_aoa_SBL, source_aoa);
        
        %% 多分辨率SBL估计
        params_MultiRes.Y = Y;
        params_MultiRes.sigma2 = noise_power;
        res_MultiRes = MultiResSBL(params_MultiRes);
        
        % 提取角度估计
        est_aoa_MultiRes = extract_angles(res_MultiRes, source_num);
        errors_MultiRes(mc) = calculate_error(est_aoa_MultiRes, source_aoa);
    end
    
    %% 计算NMSE
    NMSE_SBL(snr_idx) = mean(errors_SBL) / mean(source_aoa.^2);
    NMSE_MultiRes(snr_idx) = mean(errors_MultiRes) / mean(source_aoa.^2);
end

%% 绘制NMSE对比曲线
figure;
semilogy(SNR_range, NMSE_SBL, 'bo-', 'LineWidth', 2, 'MarkerSize', 8); hold on;
semilogy(SNR_range, NMSE_MultiRes, 'rs-', 'LineWidth', 2, 'MarkerSize', 8);
grid on;
xlabel('SNR (dB)');
ylabel('NMSE');
legend('Standard SBL', 'Multi-Resolution SBL', 'Location', 'best');
title('DOA Estimation Performance Comparison');
set(gca, 'YScale', 'log');

%% 普通SBL算法
function res = StandardSBL(params)
    Y = params.Y;
    reso = params.reso;
    grid_angles = -90:reso:90;
    A = exp(-1i*(0:size(Y,1)-1)'*pi*sind(grid_angles));
    
    % 初始化
    alpha = ones(length(grid_angles), 1);
    alpha_0 = 1/params.sigma2;
    c_sigma = 1e-4;
    d_sigma = 1e-4;
    converged = false;
    iter = 0;
    
    while ~converged && iter < params.maxiter
        iter = iter + 1;
        alpha_last = alpha;
        
        % 更新后验分布
        Sigma = inv(diag(alpha) + alpha_0*(A'*A));
        mu = alpha_0 * Sigma * A' * Y;
        
        % 更新超参数
        alpha = mean(abs(mu).^2, 2) + real(diag(Sigma));
        alpha_0 = (size(Y,2)*size(Y,1) + c_sigma - 1) / ...
                 (norm(Y - A*mu,'fro')^2 + size(Y,2)*trace(A*Sigma*A') + d_sigma);
        
        % 检查收敛
        if norm(alpha - alpha_last)/norm(alpha_last) < params.tolerance
            converged = true;
        end
    end
    
    res.mu = mu;
    res.Sigma = Sigma;
    res.grid_angles = grid_angles;
end

%% 多分辨率SBL算法
function res = MultiResSBL(params)
    Y = params.Y;
    
    %% 第一阶段: 粗网格SBL
    coarse_grid = -90:params.coarse_reso:90;
    coarse_A = exp(-1i*(0:size(Y,1)-1)'*pi*sind(coarse_grid));
    
    % 初始化
    alpha = ones(length(coarse_grid), 1);
    alpha_0 = 1/params.sigma2;
    c_sigma = 1e-4;
    d_sigma = 1e-4;
    
    for iter = 1:params.maxiter
        alpha_last = alpha;
        
        % 更新后验分布
        Sigma = inv(diag(alpha) + alpha_0*(coarse_A'*coarse_A));
        mu = alpha_0 * Sigma * coarse_A' * Y;
        
        % 更新超参数
        alpha = mean(abs(mu).^2, 2) + real(diag(Sigma));
        alpha_0 = (size(Y,2)*size(Y,1) + c_sigma - 1) / ...
                 (norm(Y - coarse_A*mu,'fro')^2 + size(Y,2)*trace(coarse_A*Sigma*coarse_A') + d_sigma);
        
        % 检查收敛
        if norm(alpha - alpha_last)/norm(alpha_last) < params.tolerance
            break;
        end
    end
    
    %% 识别需要细化的区域
    xpower = mean(abs(mu).^2,2) + real(diag(Sigma));
    xpower = xpower/max(xpower);
    [~, peak_locs] = findpeaks(xpower, 'MinPeakHeight', 0.1);
    
    % 确定精网格范围
    fine_windows = [];
    for k = 1:length(peak_locs)
        center_angle = coarse_grid(peak_locs(k));
        fine_start = max(-90, center_angle - params.window_size/2);
        fine_end = min(90, center_angle + params.window_size/2);
        fine_windows = [fine_windows; fine_start, fine_end];
    end
    
    %% 第二阶段: 精网格SBL
    % 构建混合网格
    fine_grid = [];
    for w = 1:size(fine_windows,1)
        fine_grid = [fine_grid, fine_windows(w,1):params.fine_reso:fine_windows(w,2)];
    end
    
    % 合并网格
    full_grid = unique([coarse_grid, fine_grid]);
    full_A = exp(-1i*(0:size(Y,1)-1)'*pi*sind(full_grid));
    
    % 重新初始化
    alpha = ones(length(full_grid),1);
    alpha_0 = 1/params.sigma2;
    
    for iter = 1:params.maxiter
        alpha_last = alpha;
        
        % 更新后验分布
        Sigma = inv(diag(alpha) + alpha_0*(full_A'*full_A));
        mu = alpha_0 * Sigma * full_A' * Y;
        
        % 更新超参数
        alpha = mean(abs(mu).^2,2) + real(diag(Sigma));
        alpha_0 = (size(Y,2)*size(Y,1) + c_sigma - 1) / ...
                 (norm(Y - full_A*mu,'fro')^2 + size(Y,2)*trace(full_A*Sigma*full_A') + d_sigma);
        
        % 检查收敛
        if norm(alpha - alpha_last)/norm(alpha_last) < params.tolerance
            break;
        end
    end
    
    res.mu = mu;
    res.Sigma = Sigma;
    res.grid_angles = full_grid;
end

%% 辅助函数：提取估计角度
function est_aoa = extract_angles(res, source_num)
    xpower = mean(abs(res.mu).^2, 2) + real(diag(res.Sigma));
    xpower = xpower / max(xpower);
    [~, est_peaks] = findpeaks(xpower, 'SortStr', 'descend', 'NPeaks', source_num);
    est_aoa = res.grid_angles(est_peaks);
end

%% 辅助函数：计算角度误差
function error = calculate_error(est_aoa, true_aoa)
    error = 0;
    for k = 1:length(true_aoa)
        [~, idx] = min(abs(est_aoa - true_aoa(k)));
        error = error + (est_aoa(idx) - true_aoa(k))^2;
        est_aoa(idx) = Inf; % 避免重复匹配
    end
    error = error / length(true_aoa);
end