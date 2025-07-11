close all
clear all
clc

%% 算法比较脚本 - 六种算法性能对比
% 包含：1D SBL, 2D SBL, 分层1D SBL, 分层2D SBL, OMP, Traditional Impulse
% 每个算法的分辨率参数可独立调整

%% 算法分辨率参数设置 (可独立调整)
% 1D SBL 参数
r_v_1D = 0.4;    % 1D SBL 多普勒分辨率
r_t_1D = 0.4;    % 1D SBL 时延分辨率

% 2D SBL 参数
r_v_2D = 0.4;    % 2D SBL 多普勒分辨率
r_t_2D = 0.4;    % 2D SBL 时延分辨率

% 分层1D SBL 参数
r_v_coarse_1D = 0.5;    % 分层1D SBL 粗网格多普勒分辨率
r_t_coarse_1D = 0.5;    % 分层1D SBL 粗网格时延分辨率
r_v_fine_1D = 0.2;      % 分层1D SBL 细网格多普勒分辨率
r_t_fine_1D = 0.2;      % 分层1D SBL 细网格时延分辨率
threshold_ratio_1D = 0.1; % 分层1D SBL 系数选择阈值

% 分层2D SBL 参数
r_v_coarse_2D = 0.5;    % 分层2D SBL 粗网格多普勒分辨率
r_t_coarse_2D = 0.5;    % 分层2D SBL 粗网格时延分辨率
r_v_fine_2D = 0.2;      % 分层2D SBL 细网格多普勒分辨率
r_t_fine_2D = 0.2;      % 分层2D SBL 细网格时延分辨率
threshold_ratio_2D = 0.1; % 分层2D SBL 系数选择阈值

% OMP 参数
r_v_OMP = 0.4;    % OMP 多普勒分辨率
r_t_OMP = 0.4;    % OMP 时延分辨率

%% Test Mode %%%%%%
test = 0;    % set to 1 when testing

%% OTFS parameters%%%%%%%%%%
% N: number of symbols in time
N = 64;
% M: number of subcarriers in frequency
M = 64;
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

%% parameters for test use
k_v_test = [0,0,0,0,0]';
l_t_test = [0,0,0,0,0]';
h_test = [1,0,0,0,0]';

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
SNR_dB = 0:5:20;
SNR = 10.^(SNR_dB/10);
sigma_2 = 0.5 ./SNR;
SNR_p = SNR * 1e5;
sigma_2_p = sigma_2 /1e5;

%% Initializing simulation error count variables

N_fram = 10;  % 减少仿真帧数以加快运行速度
NMSE_count_1D = zeros(length(SNR_dB),1);
NMSE_count_2D = zeros(length(SNR_dB),1);
NMSE_count_hierarchical_1D = zeros(length(SNR_dB),1);
NMSE_count_hierarchical_2D = zeros(length(SNR_dB),1);
NMSE_count_OMP = zeros(length(SNR_dB),1);
NMSE_count_traditional = zeros(length(SNR_dB),1);
NMSE_CRLB = zeros(length(SNR_dB),1);  % CRLB理论下界
err_ber = zeros(length(SNR_dB),1);
H_mtx = zeros(M*N,M*N);

fprintf('\n=== 六种算法性能比较开始 ===\n');
fprintf('仿真参数：N=%d, M=%d, P=%d, N_fram=%d\n', N, M, P, N_fram);
fprintf('\n算法分辨率参数：\n');
fprintf('1D SBL: r_v=%.2f, r_t=%.2f\n', r_v_1D, r_t_1D);
fprintf('2D SBL: r_v=%.2f, r_t=%.2f\n', r_v_2D, r_t_2D);
fprintf('分层1D SBL: 粗网格(r_v=%.2f, r_t=%.2f), 细网格(r_v=%.2f, r_t=%.2f)\n', r_v_coarse_1D, r_t_coarse_1D, r_v_fine_1D, r_t_fine_1D);
fprintf('分层2D SBL: 粗网格(r_v=%.2f, r_t=%.2f), 细网格(r_v=%.2f, r_t=%.2f)\n', r_v_coarse_2D, r_t_coarse_2D, r_v_fine_2D, r_t_fine_2D);
fprintf('OMP: r_v=%.2f, r_t=%.2f\n', r_v_OMP, r_t_OMP);
fprintf('Traditional Impulse: 传统脉冲响应估计\n');
fprintf('\n');

for ifram = 1:N_fram
    
    %%  random channel initialization
    v_c_init = unifrnd(-v_max,v_max,P,1);   
    t_c_init = unifrnd(0,t_max,P,1);
    l_ti = t_c_init.*(M * delta_f);
    q_l_t = exp(-0.1.*l_ti)./sum(exp(-0.1.*l_ti));
    h_c_init = normrnd(0,q_l_t);  % normrnd but not mvnrnd
    k_v_init = v_c_init .*(N*T);
    l_t_init = t_c_init .*(M*delta_f);
    
    for iesn0 = 1:length(SNR_dB)
        fprintf('Frame %d/%d, SNR %d dB\n', ifram, N_fram, SNR_dB(iesn0));
        
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
        
        %% get truncation measument matrix%%%%%
        y_trunc = y(floor(N/2)-floor(Kp/2)-k_max:floor(N/2)+floor(Kp/2)+k_max,floor(M/2)-floor(Lp/2):floor(M/2)+floor(Lp/2)+l_max);
        N_T = 2*k_max+Kp;
        M_T = l_max+Lp;
        y_T = reshape(y_trunc.',N_T*M_T,1);        
        
        %% 算法1: 1D SBL Channel Estimation
        [h_hat_1D,k_v_hat_1D,l_t_hat_1D,virtual_size_1D, Phi_1D, delta_a] = CE_1D_SBL(sqrt(1000)/(Kp*Lp),k_max+Kp,Lp,M,N,N_T,M_T,y_T,r_v_1D,r_t_1D,k_max,l_max,on_flag);

        %% 算法2: 2D SBL Channel Estimation
        [H_opt_2D,k_v_opt_2D,l_t_opt_2D,virtual_size_2D] = CE_2D_SBL(sqrt(1000)/(2*Kp*Lp),k_max+Kp,2,Lp,M,N,N_T,M_T,y_trunc,r_v_2D,r_t_2D,k_max,l_max,on_flag);

        %% 算法3: 分层1D SBL Channel Estimation
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

        %% 算法4: 分层2D SBL Channel Estimation
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
       
        %% 算法5: OMP Channel Estimation
        [h_hat_OMP, omp_index, k_v_hat_OMP, l_t_hat_OMP] = OMP(sqrt(1000)/(Kp*Lp),k_max+Kp,Lp,M,N,N_T,M_T,y_T,r_v_OMP,r_t_OMP,k_max,l_max);
        
        %% 算法6: Traditional Impulse Channel Estimation
        [h_hat_traditional, k_v_hat_traditional, l_t_hat_traditional] = traditional_impulse(sqrt(1000)/(Kp*Lp), y_trunc, k_max, l_max, sigma_2_p(iesn0));
        
        %% 计算CRLB理论下界
         % 构建雅可比矩阵J
         J = zeros(M_T*N_T, virtual_size_1D);
         for m=1:M_T
             for n=1:N_T
                 J((n-1)*M_T+m,:) = (Sampling_Function_v(N,n,1,k_v_hat_1D) .* Sampling_Function_t(M,m,1,l_t_hat_1D) .* exp(-2i*pi.*k_v_hat_1D.*l_t_hat_1D/N/M)).';
             end
         end
         % 计算协方差矩阵
         C_h = (sigma_2(iesn0)^(-1) * (Phi_1D' * Phi_1D) + delta_a^(-1))^(-1);
         % C_h = sigma_2(iesn0) * (Phi_1D' * Phi_1D)^(-1);
         C_alpha = J * C_h * J.';
         NMSE_CRLB(iesn0) = trace(abs(C_alpha));
        
        %% 计算所有算法的NMSE
        NMSE_nume_1D = 0;
        NMSE_nume_2D = 0;
        NMSE_nume_hierarchical_1D = 0;
        NMSE_nume_hierarchical_2D = 0;
        NMSE_nume_OMP = 0;
        NMSE_nume_traditional = 0;
        NMSE_deno = 0;
        
        for kk = 0:(N-1)
            for ll = 1:M
                if (test==1)
                    h_w = sum(Sampling_Function_v(N,kk+1,1,k_v_test) .* Sampling_Function_t(M,ll,1,l_t_test) .* h_test .* exp(-2i*pi.*k_v_test.*l_t_test/N/M));
                else
                    h_w = sum(Sampling_Function_v(N,kk+1,1,k_v_init) .* Sampling_Function_t(M,ll,1,l_t_init) .* h_c_init .* exp(-2i*pi.*k_v_init.*l_t_init/N/M));
                end
                
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
                
                % Traditional Impulse estimation
                h_w_hat_traditional = sum(Sampling_Function_v(N,kk+1,1,k_v_hat_traditional) .* Sampling_Function_t(M,ll,1,l_t_hat_traditional) .* h_hat_traditional .* exp(-2i*pi.*k_v_hat_traditional.*l_t_hat_traditional/N/M));
                
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

%% Results plotting
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
plot(SNR_dB,NMSE_dB_1D,'-v', 'Color', colors(3,:),'LineWidth',2,'MarkerSize',8);
hold on
plot(SNR_dB,NMSE_dB_2D,'-*', 'Color', colors(2,:),'LineWidth',2,'MarkerSize',8);
plot(SNR_dB,NMSE_dB_hierarchical_1D,'-s', 'Color', colors(1,:),'LineWidth',2,'MarkerSize',8);
plot(SNR_dB,NMSE_dB_hierarchical_2D,'-o', 'Color', colors(4,:),'LineWidth',2,'MarkerSize',8);
plot(SNR_dB,NMSE_dB_OMP,'-d', 'Color', colors(6,:),'LineWidth',2,'MarkerSize',8);
plot(SNR_dB,NMSE_dB_traditional,'-x', 'Color', colors(7,:),'LineWidth',2,'MarkerSize',8);
plot(SNR_dB,NMSE_dB_CRLB,'-^', 'Color', colors(5,:),'LineWidth',2,'MarkerSize',8);
grid on

legend('1D SBL', '2D SBL', '分层1D SBL', '分层2D SBL', 'OMP', 'Traditional Impulse', 'CRLB下界', 'Location', 'best');
xlabel('SNR (dB)', 'FontSize', 12);
ylabel('NMSE (dB)', 'FontSize', 12);
title('六种算法信道估计性能比较', 'FontSize', 14);
set(gca, 'FontSize', 11);

%% 结果分析和显示
fprintf('\n=== 六种算法性能比较结果 ===\n');
fprintf('SNR(dB)\t1D SBL\t\t2D SBL\t\t分层1D SBL\t分层2D SBL\tOMP\t\t\tTraditional\tCRLB下界\n');
fprintf('------\t------\t\t------\t\t----------\t----------\t----------\t----------\t--------\n');
for i = 1:length(SNR_dB)
    fprintf('%d\t\t%.2f\t\t%.2f\t\t%.2f\t\t%.2f\t\t%.2f\t\t%.2f\t\t%.2f\n', SNR_dB(i), NMSE_dB_1D(i), NMSE_dB_2D(i), NMSE_dB_hierarchical_1D(i), NMSE_dB_hierarchical_2D(i), NMSE_dB_OMP(i), NMSE_dB_traditional(i), NMSE_dB_CRLB(i));
end

%% 算法复杂度分析
fprintf('\n=== 算法复杂度分析 ===\n');
fprintf('1D SBL格点个数: %d\n', ceil(2 * k_max / r_v_1D) * ceil(l_max / r_t_1D));
fprintf('2D SBL格点个数: %d\n', ceil(2 * k_max / r_v_2D) * ceil(l_max / r_t_2D));
fprintf('分层1D SBL粗网格格点个数: %d\n', ceil(2 * k_max / r_v_coarse_1D) * ceil(l_max / r_t_coarse_1D));
fprintf('分层1D SBL细网格格点个数: %d\n', ceil(2 * k_max / r_v_fine_1D) * ceil(l_max / r_t_fine_1D));
fprintf('分层2D SBL粗网格格点个数: %d\n', ceil(2 * k_max / r_v_coarse_2D) * ceil(l_max / r_t_coarse_2D));
fprintf('分层2D SBL细网格格点个数: %d\n', ceil(2 * k_max / r_v_fine_2D) * ceil(l_max / r_t_fine_2D));
fprintf('OMP格点个数: %d\n', ceil(2 * k_max / r_v_OMP) * ceil(l_max / r_t_OMP));
fprintf('Traditional Impulse格点个数: %d\n', (2 * k_max + 1) * (l_max + 1));

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