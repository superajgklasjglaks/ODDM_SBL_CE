% Massive MIMO 1D SBL 方向角估计算法
% 简化版本，只包含1D SBL算法
% 作者: AI Assistant
% 日期: 2024

clc; clear; close all;

%% 系统参数设置
% Massive MIMO系统参数
Nt = 32;                    % 发射天线数
Nr = 32;                    % 接收天线数
P = 3;                      % 路径数
L_pilot = 16;               % 导频长度

% 角度参数
aod_range = [-pi/3, pi/3];  % 发射角范围 (弧度)
aoa_range = [-pi/3, pi/3];  % 接收角范围 (弧度)

% 算法参数
resolution_1D = 0.1;       % 1D网格分辨率
virtual_factor = 4;         % 虚拟网格因子

% 仿真参数
num_frames = 100;           % 仿真帧数
SNR_dB = 10;                % 信噪比 (dB)

%% 生成真实角度
% 随机生成发射角和接收角
aod_true = aod_range(1) + (aod_range(2) - aod_range(1)) * rand(P, 1);
aoa_true = aoa_range(1) + (aoa_range(2) - aoa_range(1)) * rand(P, 1);

% 生成真实信道系数
h_true = (randn(P, 1) + 1j * randn(P, 1)) / sqrt(2);

fprintf('真实发射角 (度): ');
fprintf('%.2f ', aod_true * 180 / pi);
fprintf('\n');
fprintf('真实接收角 (度): ');
fprintf('%.2f ', aoa_true * 180 / pi);
fprintf('\n\n');

%% 1D SBL算法参数
% 生成角度网格
aod_grid = aod_range(1):resolution_1D:aod_range(2);
aoa_grid = aoa_range(1):resolution_1D:aoa_range(2);
G_aod = length(aod_grid);
G_aoa = length(aoa_grid);

fprintf('发射角网格点数: %d\n', G_aod);
fprintf('接收角网格点数: %d\n', G_aoa);

%% 初始化结果存储
NMSE_aod_1D = zeros(num_frames, 1);
NMSE_aoa_1D = zeros(num_frames, 1);
NMSE_h_1D = zeros(num_frames, 1);

%% Monte Carlo仿真
fprintf('\n开始1D SBL算法仿真...\n');
tic;

for frame = 1:num_frames
    if mod(frame, 20) == 0
        fprintf('完成帧数: %d/%d\n', frame, num_frames);
    end
    
    %% 生成导频信号和接收信号
    % 导频矩阵 (随机导频)
    X_pilot = (randn(Nt, L_pilot) + 1j * randn(Nt, L_pilot)) / sqrt(2);
    
    % 生成信道矩阵
    H = zeros(Nr, Nt);
    for p = 1:P
        % 发射导向向量
        a_tx = exp(1j * (0:Nt-1)' * pi * sin(aod_true(p)));
        % 接收导向向量
        a_rx = exp(1j * (0:Nr-1)' * pi * sin(aoa_true(p)));
        % 信道贡献
        H = H + h_true(p) * a_rx * a_tx';
    end
    
    % 接收信号
    Y_pilot = H * X_pilot;
    
    % 添加噪声
    noise_power = 10^(-SNR_dB/10);
    noise = sqrt(noise_power/2) * (randn(Nr, L_pilot) + 1j * randn(Nr, L_pilot));
    Y_pilot_noisy = Y_pilot + noise;
    
    % 向量化接收信号
    y = Y_pilot_noisy(:);
    
    %% 1D SBL算法
    try
        [h_hat_1D, aod_hat_1D, aoa_hat_1D] = MIMO_CE_1D_SBL(...
            y, X_pilot, Nt, Nr, L_pilot, P, ...
            aod_grid, aoa_grid, virtual_factor);
        
        % 计算NMSE
        if ~isempty(aod_hat_1D) && ~isempty(aoa_hat_1D)
            % 角度匹配 (找到最接近的估计值)
            aod_error = min(abs(aod_hat_1D - aod_true'));
            aoa_error = min(abs(aoa_hat_1D - aoa_true'));
            
            NMSE_aod_1D(frame) = mean(aod_error.^2) / mean(aod_true.^2);
            NMSE_aoa_1D(frame) = mean(aoa_error.^2) / mean(aoa_true.^2);
            
            if ~isempty(h_hat_1D)
                h_error = min(abs(h_hat_1D - h_true'));
                NMSE_h_1D(frame) = mean(h_error.^2) / mean(abs(h_true).^2);
            else
                NMSE_h_1D(frame) = 1;
            end
        else
            NMSE_aod_1D(frame) = 1;
            NMSE_aoa_1D(frame) = 1;
            NMSE_h_1D(frame) = 1;
        end
        
    catch ME
        fprintf('第%d帧1D SBL算法出错: %s\n', frame, ME.message);
        NMSE_aod_1D(frame) = 1;
        NMSE_aoa_1D(frame) = 1;
        NMSE_h_1D(frame) = 1;
    end
end

simulation_time = toc;
fprintf('仿真完成，总耗时: %.2f秒\n\n', simulation_time);

%% 结果统计
fprintf('=== 1D SBL算法性能统计 ===\n');
fprintf('发射角NMSE (dB): %.2f\n', 10*log10(mean(NMSE_aod_1D)));
fprintf('接收角NMSE (dB): %.2f\n', 10*log10(mean(NMSE_aoa_1D)));
fprintf('信道系数NMSE (dB): %.2f\n', 10*log10(mean(NMSE_h_1D)));
fprintf('平均每帧耗时: %.4f秒\n', simulation_time/num_frames);

%% 绘制结果
figure('Position', [100, 100, 1200, 400]);

% 发射角NMSE
subplot(1, 3, 1);
plot(1:num_frames, 10*log10(NMSE_aod_1D), 'b-', 'LineWidth', 1.5);
hold on;
yline(10*log10(mean(NMSE_aod_1D)), 'r--', 'LineWidth', 2);
xlabel('仿真帧数');
ylabel('发射角NMSE (dB)');
title('1D SBL - 发射角估计性能');
grid on;
legend('瞬时NMSE', sprintf('平均NMSE: %.2f dB', 10*log10(mean(NMSE_aod_1D))), 'Location', 'best');

% 接收角NMSE
subplot(1, 3, 2);
plot(1:num_frames, 10*log10(NMSE_aoa_1D), 'g-', 'LineWidth', 1.5);
hold on;
yline(10*log10(mean(NMSE_aoa_1D)), 'r--', 'LineWidth', 2);
xlabel('仿真帧数');
ylabel('接收角NMSE (dB)');
title('1D SBL - 接收角估计性能');
grid on;
legend('瞬时NMSE', sprintf('平均NMSE: %.2f dB', 10*log10(mean(NMSE_aoa_1D))), 'Location', 'best');

% 信道系数NMSE
subplot(1, 3, 3);
plot(1:num_frames, 10*log10(NMSE_h_1D), 'm-', 'LineWidth', 1.5);
hold on;
yline(10*log10(mean(NMSE_h_1D)), 'r--', 'LineWidth', 2);
xlabel('仿真帧数');
ylabel('信道系数NMSE (dB)');
title('1D SBL - 信道系数估计性能');
grid on;
legend('瞬时NMSE', sprintf('平均NMSE: %.2f dB', 10*log10(mean(NMSE_h_1D))), 'Location', 'best');

sgtitle('Massive MIMO 1D SBL 方向角估计性能', 'FontSize', 14, 'FontWeight', 'bold');

%% 显示最终估计结果
fprintf('\n=== 最后一帧的估计结果 ===\n');
if exist('aod_hat_1D', 'var') && ~isempty(aod_hat_1D)
    fprintf('估计发射角 (度): ');
    fprintf('%.2f ', aod_hat_1D * 180 / pi);
    fprintf('\n');
end
if exist('aoa_hat_1D', 'var') && ~isempty(aoa_hat_1D)
    fprintf('估计接收角 (度): ');
    fprintf('%.2f ', aoa_hat_1D * 180 / pi);
    fprintf('\n');
end
if exist('h_hat_1D', 'var') && ~isempty(h_hat_1D)
    fprintf('估计信道系数幅度: ');
    fprintf('%.4f ', abs(h_hat_1D));
    fprintf('\n');
end

fprintf('\n1D SBL算法仿真完成！\n');