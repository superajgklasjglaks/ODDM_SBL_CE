close all
clear all
clc
%rng(128)

%% 测试模式 %%%%%%
test = 0;    % 测试时设置为1

%% OTFS参数%%%%%%%%%%
% N: 时域符号数
N = 32;
% M: 频域子载波数
M = 32;
% M_mod: QAM星座图大小
M_mod = 4;
M_bits = log2(M_mod);
% 每个数据符号的平均能量
eng_sqrt = (M_mod==2)+(M_mod~=2)*sqrt((M_mod-1)/6*(2^2));

%% 时延-多普勒网格符号放置

% 时间和频率资源
car_fre = 28*10^9;  % 载波频率
delta_f = 625*10^3; % 子载波间隔: 15 KHz
T = 1/delta_f;     % OTFS帧中一个时间符号的持续时间
k_max = 4;         % 系统设置
l_max = 4;
t_max = l_max / M / delta_f;
v_max = k_max / N / T;

%% 分层SBL参数
% 粗网格参数
r_v_coarse = 0.5;    % 多普勒粗分辨率
r_t_coarse = 0.5;    % 时延粗分辨率
% 细网格参数
r_v_fine = 0.2;      % 多普勒细分辨率
r_t_fine = 0.2;      % 时延细分辨率
% 网格细化的选择阈值
threshold_ratio = 0.1; % 选择前10%的估计系数
P = 4;
on_flag = 0;      % 在网格上时设置为0

%% 测试用参数
k_v_test = [0,0,0,0,0]';
l_t_test = [0,0,0,0,0]';
h_test = [1,0,0,0,0]';

%% 导频符号放置
Kp = 1;
Lp = 1;
x_kp = floor(N/2);
x_lp = floor(M/2);
% OTFS时延-多普勒域数据符号在2D网格中的数据位置
data_grid = ones(N,M);
data_grid(x_kp-floor(Kp/2)-2*k_max:x_kp+floor(Kp/2)+2*k_max,x_lp-floor(Lp/2)-l_max:x_lp+floor(Lp/2)+l_max)=0;
% 每帧符号数
N_syms_perfram = sum(sum(data_grid));

% 每帧比特数
N_bits_perfram = N_syms_perfram*M_bits;

 
% 信噪比和噪声方差
% SNR = P/\sigma^2; P: 传输字母表的平均功率
SNR_dB = 0:5:20;
SNR = 10.^(SNR_dB/10);
sigma_2 = 0.5 ./SNR;
SNR_p = SNR * 1e5;
sigma_2_p = sigma_2 /1e5;

%% 初始化仿真误差计数变量

N_fram = 10;
NMSE_count_coarse = zeros(length(SNR_dB),1);
NMSE_count_hierarchical = zeros(length(SNR_dB),1);
err_ber = zeros(length(SNR_dB),1);
H_mtx = zeros(M*N,M*N);

for ifram = 1:N_fram
    %%  随机信道初始化
    v_c_init = unifrnd(-v_max,v_max,P,1);   
    t_c_init = unifrnd(0,t_max,P,1);
    l_ti = t_c_init.*(M * delta_f);
    q_l_t = exp(-0.1.*l_ti)./sum(exp(-0.1.*l_ti));
    h_c_init = normrnd(0,q_l_t);  % 使用normrnd而不是mvnrnd
    k_v_init = v_c_init .*(N*T);
    l_t_init = t_c_init .*(M*delta_f);
    
    for iesn0 = 1:length(SNR_dB)
        pow_prof = (1/P) * (ones(1,P));
        chan_coef = zeros(1,P);
        
        %% 随机输入比特生成%%%%%
        trans_info_bit = randi([0,1],N_syms_perfram*M_bits,1);
        %%2D QAM符号生成 %%%%%%%%
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
        
        %% 信道输出%%%%%
        noise_gen_re = sigma_2(iesn0) * randn(M*N,1);
        noise_gen_im = sigma_2(iesn0) * randn(M*N,1);
        noise_gen = noise_gen_re + 1i.* noise_gen_im;
        h_test_hat = h_test .* exp(-2i * pi/M/N * (k_v_test.*l_t_test));
        h_c_init_hat = h_c_init .* exp(-2i * pi/M/N * (k_v_init.*l_t_init));
        if (test==1)
            r = phi_sys * h_test;
        else
            r = phi_sys * h_c_init + noise_gen;
        end
        y = reshape(r,M,N).';
        
        %% 获取截断测量矩阵%%%%%
        y_trunc = y(floor(N/2)-floor(Kp/2)-k_max:floor(N/2)+floor(Kp/2)+k_max,floor(M/2)-floor(Lp/2):floor(M/2)+floor(Lp/2)+l_max);
        N_T = 2*k_max+Kp;
        M_T = l_max+Lp;
        y_T = reshape(y_trunc.',N_T*M_T,1);        
        
        %% 步骤1: 粗网格SBL估计
        fprintf('Frame %d, SNR %d dB: Starting coarse grid SBL...\n', ifram, SNR_dB(iesn0));
        [h_hat_coarse, k_v_hat_coarse, l_t_hat_coarse, virtual_size_coarse, Phi_coarse, delta_coarse] = ...
            CE_1D_SBL(sqrt(1000)/(Kp*Lp), k_max+Kp, Lp, M, N, N_T, M_T, y_T, r_v_coarse, r_t_coarse, k_max, l_max, on_flag);
        
        %% 步骤2: 识别重要系数及其位置
        [sorted_coeff, sorted_idx] = sort(abs(h_hat_coarse), 'descend');
        num_significant = max(1, floor(threshold_ratio * length(h_hat_coarse)));
        significant_idx = sorted_idx(1:num_significant);
        
        % 提取重要系数对应的k_v和l_t值
        k_v_significant = k_v_hat_coarse(significant_idx);
        l_t_significant = l_t_hat_coarse(significant_idx);
        
        fprintf('  Selected %d significant coefficients for refinement\n', num_significant);
        
        %% 步骤3: 在重要位置周围创建细化网格
        % 定义每个重要位置周围的细化网格边界
        refined_k_v = [];
        refined_l_t = [];
        refined_regions = [];
        
        for i = 1:length(k_v_significant)
            % 定义重要系数周围的局部区域
            k_center = k_v_significant(i);
            l_center = l_t_significant(i);
            
            % 在此中心周围创建细网格
            k_range = k_center + (-r_v_coarse/2:r_v_fine:r_v_coarse/2);
            l_range = l_center + (-r_t_coarse/2:r_t_fine:r_t_coarse/2);
            
            % 确保边界在有效范围内
            k_range = k_range(k_range >= -k_max & k_range <= k_max);
            l_range = l_range(l_range >= 0 & l_range <= l_max);
            
            % 为此区域创建网格
            [K_mesh, L_mesh] = meshgrid(k_range, l_range);
            region_k = K_mesh(:);
            region_l = L_mesh(:);
            
            refined_k_v = [refined_k_v; region_k];
            refined_l_t = [refined_l_t; region_l];
            refined_regions = [refined_regions; i*ones(length(region_k), 1)];
        end
        
        % 去除重复项
        [refined_grid, unique_idx] = unique([refined_k_v, refined_l_t], 'rows');
        refined_k_v = refined_grid(:, 1);
        refined_l_t = refined_grid(:, 2);
        
        fprintf('  Created refined grid with %d points\n', length(refined_k_v));
        
        %% 步骤4: 在细化区域上执行细网格SBL
        if length(refined_k_v) > 0
            % 为细化网格创建自定义测量矩阵
            [h_hat_fine, k_v_hat_fine, l_t_hat_fine] = ...
                hierarchical_SBL_refined_1D(sqrt(1000)/(Kp*Lp), k_max+Kp, Lp, M, N, N_T, M_T, y_T, ...
                                        refined_k_v, refined_l_t, on_flag);
            
            fprintf('  Completed fine grid SBL\n');
        else
            % 回退到粗网格结果
            h_hat_fine = h_hat_coarse;
            k_v_hat_fine = k_v_hat_coarse;
            l_t_hat_fine = l_t_hat_coarse;
        end
        
        %% 步骤5: 计算粗网格和分层方法的NMSE
        NMSE_nume_coarse = 0;
        NMSE_nume_hierarchical = 0;
        NMSE_deno = 0;
        
        for kk = 0:(N-1)
            for ll = 1:M
                if (test==1)
                    h_w = sum(Sampling_Function_v(N,kk+1,1,k_v_test) .* Sampling_Function_t(M,ll,1,l_t_test) .* h_test .* exp(-2i*pi.*k_v_test.*l_t_test/N/M));
                else
                    h_w = sum(Sampling_Function_v(N,kk+1,1,k_v_init) .* Sampling_Function_t(M,ll,1,l_t_init) .* h_c_init .* exp(-2i*pi.*k_v_init.*l_t_init/N/M));
                end
                
                % 粗网格估计
                h_w_hat_coarse = sum(Sampling_Function_v(N,kk+1,1,k_v_hat_coarse) .* Sampling_Function_t(M,ll,1,l_t_hat_coarse) .* h_hat_coarse .* exp(-2i*pi.*k_v_hat_coarse.*l_t_hat_coarse/N/M));
                
                % 分层估计
                h_w_hat_hierarchical = sum(Sampling_Function_v(N,kk+1,1,k_v_hat_fine) .* Sampling_Function_t(M,ll,1,l_t_hat_fine) .* h_hat_fine .* exp(-2i*pi.*k_v_hat_fine.*l_t_hat_fine/N/M));
                
                NMSE_nume_coarse = NMSE_nume_coarse + abs(h_w - h_w_hat_coarse).^2;
                NMSE_nume_hierarchical = NMSE_nume_hierarchical + abs(h_w - h_w_hat_hierarchical).^2;
                NMSE_deno = NMSE_deno + abs(h_w)^2;
            end 
        end
        
        % 不需要CRLB计算
        
        NMSE_count_coarse(iesn0) = NMSE_count_coarse(iesn0) + NMSE_nume_coarse / (NMSE_deno * N_fram);
        NMSE_count_hierarchical(iesn0) = NMSE_count_hierarchical(iesn0) + NMSE_nume_hierarchical / (NMSE_deno * N_fram);
        
        fprintf('  Coarse NMSE: %.4f, Hierarchical NMSE: %.4f\n', NMSE_nume_coarse/NMSE_deno, NMSE_nume_hierarchical/NMSE_deno);
        
        display(ifram,'ifram');
        display(iesn0, 'iesn0');
    end
end

%%结果绘图
NMSE_dB_coarse = 10 * log10(NMSE_count_coarse);
NMSE_dB_hierarchical = 10 * log10(NMSE_count_hierarchical);

colors = [0,107,182; %?1
          118,41,133; %?2
          234,174,31; %?3
          215,94,59; %?4
          184,125,182; %??5
          71,90,40; %?6
          161,27,30]; %?7
colors = colors/256;

figure()
plot(SNR_dB,NMSE_dB_coarse,'-s', 'Color', colors(1,:),'LineWidth',1.5,'MarkerSize',8);
hold on
plot(SNR_dB,NMSE_dB_hierarchical,'-v', 'Color', colors(3,:),'LineWidth',1.5,'MarkerSize',8);
grid on

legend('Coarse Grid SBL', 'Hierarchical SBL');
xlabel('SNR(dB)')
ylabel('NMSE(dB)')
title('Hierarchical Sparse Bayesian Learning Channel Estimation Performance')

fprintf('\n=== Hierarchical SBL Channel Estimation Results ===\n');
fprintf('SNR(dB)\tCoarse NMSE(dB)\tHierarchical NMSE(dB)\tImprovement(dB)\n');
for i = 1:length(SNR_dB)
    improvement = NMSE_dB_coarse(i) - NMSE_dB_hierarchical(i);
    fprintf('%d\t\t%.2f\t\t\t%.2f\t\t\t\t%.2f\n', SNR_dB(i), NMSE_dB_coarse(i), NMSE_dB_hierarchical(i), improvement);
end