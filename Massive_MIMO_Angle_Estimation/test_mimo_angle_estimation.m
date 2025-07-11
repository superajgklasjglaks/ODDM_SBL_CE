% 测试Massive MIMO角度估计代码
close all
clear all
clc

fprintf('开始测试Massive MIMO角度估计代码...\n');

try
    % 设置简化的测试参数
    Nt = 8;  % 减少天线数以加快测试
    Nr = 8;
    P = 2;   % 减少路径数
    N_fram = 2;  % 减少仿真帧数
    
    % 角度参数
    aod_max = pi/6;  % 30度
    aoa_max = pi/6;
    
    % 分辨率参数
    r_aod = 0.05;
    r_aoa = 0.05;
    
    % SNR设置
    SNR_dB = [0, 10];
    SNR = 10.^(SNR_dB/10);
    sigma_2 = 0.5 ./ SNR;
    
    fprintf('测试参数设置完成\n');
    
    % 测试阵列响应函数
    fprintf('测试阵列响应函数...\n');
    aod_test = 0.1;
    aoa_test = 0.2;
    
    w_tx = MIMO_Array_Response_Tx(Nt, 1, 1, aod_test);
    w_rx = MIMO_Array_Response_Rx(Nr, 1, 1, aoa_test);
    
    fprintf('发射阵列响应: %.4f + %.4fi\n', real(w_tx), imag(w_tx));
    fprintf('接收阵列响应: %.4f + %.4fi\n', real(w_rx), imag(w_rx));
    
    % 测试角度网格生成
    fprintf('测试角度网格生成...\n');
    N_aod = ceil(2 * aod_max / r_aod);
    M_aoa = ceil(2 * aoa_max / r_aoa);
    
    [aod_bar, aoa_bar] = MIMO_First_Order_Linear_Approximation(N_aod, M_aoa, aod_max, r_aod, r_aoa);
    fprintf('生成角度网格: %d x %d = %d 个点\n', N_aod, M_aoa, length(aod_bar));
    
    % 测试简单的1D SBL
    fprintf('测试1D SBL算法...\n');
    
    % 生成测试数据
    aod_true = [0.1, -0.15];
    aoa_true = [0.2, -0.1];
    h_true = [1, 0.5];
    
    % 构建测量矩阵
    M_T = 4;
    N_T = 4;
    y_T = zeros(M_T * N_T, 1);
    
    for p = 1:length(aod_true)
        for nnn = 0:(N_T-1)
            for mmm = 1:M_T
                matrix_idx = nnn * M_T + mmm;
                y_T(matrix_idx) = y_T(matrix_idx) + h_true(p) * ...
                    MIMO_Array_Response_Tx(Nt, nnn, 3, aod_true(p)) * ...
                    MIMO_Array_Response_Rx(Nr, mmm-1, 1, aoa_true(p));
            end
        end
    end
    
    % 添加噪声
    noise = (sigma_2(1) * randn(size(y_T)) + 1j * sigma_2(1) * randn(size(y_T)));
    y_T = y_T + noise;
    
    % 调用1D SBL
    [h_hat, aod_hat, aoa_hat, virtual_size, phi_trunc, delta_a] = ...
        MIMO_CE_1D_SBL(1, 4, 2, Nr, Nt, N_T, M_T, y_T, r_aod, r_aoa, aod_max, aoa_max, 0);
    
    fprintf('1D SBL完成，估计了 %d 个系数\n', length(h_hat));
    fprintf('最大系数幅度: %.4f\n', max(abs(h_hat)));
    
    % 测试OMP
    fprintf('测试OMP算法...\n');
    [h_hat_omp, omp_index, aod_hat_omp, aoa_hat_omp] = ...
        MIMO_OMP(1, 4, 2, Nr, Nt, N_T, M_T, y_T, r_aod, r_aoa, aod_max, aoa_max);
    
    fprintf('OMP完成，选择了 %d 个原子\n', length(omp_index));
    
    % 测试传统波束成形
    fprintf('测试传统波束成形...\n');
    y_trunc = reshape(y_T, M_T, N_T).';
    [h_hat_trad, aod_hat_trad, aoa_hat_trad] = ...
        MIMO_traditional_beamforming(1, y_trunc, 2, 2, sigma_2(1));
    
    fprintf('传统波束成形完成，估计了 %d 个系数\n', length(h_hat_trad));
    
    fprintf('\n=== 所有测试通过! ===\n');
    fprintf('代码可以正常运行，可以执行完整的仿真。\n');
    
catch ME
    fprintf('\n=== 测试失败 ===\n');
    fprintf('错误信息: %s\n', ME.message);
    fprintf('错误位置: %s (第 %d 行)\n', ME.stack(1).name, ME.stack(1).line);
    
    % 显示详细的错误堆栈
    for i = 1:length(ME.stack)
        fprintf('  %d. %s (第 %d 行)\n', i, ME.stack(i).name, ME.stack(i).line);
    end
end