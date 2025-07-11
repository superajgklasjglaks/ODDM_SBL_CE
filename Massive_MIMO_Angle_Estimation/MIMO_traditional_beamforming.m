function [h_hat, aod_hat, aoa_hat] = MIMO_traditional_beamforming(x_p, y_trunc, aod_range, aoa_range, sigma_2_p)
% Massive MIMO 传统波束成形角度估计函数
% 模仿ODDM中的traditional_impulse
% 输入:
%   x_p: 导频功率
%   y_trunc: 截断接收信号矩阵
%   aod_range, aoa_range: 角度搜索范围
%   sigma_2_p: 噪声功率
% 输出:
%   h_hat: 估计的信道系数
%   aod_hat, aoa_hat: 估计的角度

[M_T, N_T] = size(y_trunc);

% 简单的网格搜索
aod_grid = linspace(-pi/3, pi/3, 2*aod_range+1);
aoa_grid = linspace(-pi/3, pi/3, 2*aoa_range+1);

virtual_size = length(aod_grid) * length(aoa_grid);
h_hat = zeros(virtual_size, 1);
aod_hat = zeros(virtual_size, 1);
aoa_hat = zeros(virtual_size, 1);

% 天线间距
d_lambda = 0.5;

index = 1;
for i = 1:length(aod_grid)
    for j = 1:length(aoa_grid)
        aod = aod_grid(i);
        aoa = aoa_grid(j);
        
        % 构建导向向量
        at = exp(1j * 2 * pi * d_lambda * (0:N_T-1)' * sin(aod)) / sqrt(N_T);
        ar = exp(1j * 2 * pi * d_lambda * (0:M_T-1)' * sin(aoa)) / sqrt(M_T);
        
        % 波束成形
        beamformed_signal = ar' * y_trunc * at;
        
        % 估计信道系数
        h_hat(index) = beamformed_signal / x_p;
        aod_hat(index) = aod;
        aoa_hat(index) = aoa;
        
        index = index + 1;
    end
end

% 添加噪声抑制
noise_threshold = sqrt(sigma_2_p) * 3;
h_hat(abs(h_hat) < noise_threshold) = 0;

end