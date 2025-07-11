close all
clear all
clc

%% �㷨������·������P�仯����
% �����꣺·������P = 3:3:15 (��3, 6, 9, 12, 15)
% �����꣺NMSE����
% ������1D SBL, 2D SBL, �ֲ�1D SBL, �ֲ�2D SBL, OMP, Traditional Impulse

%% ·��������������
P_values = 3:3:15;  % ·��������3, 6, 9, 12, 15
num_P_values = length(P_values);

%% �㷨�ֱ��ʲ������� (�̶�)
% 1D SBL ����
r_v_1D = 0.4;    % 1D SBL �����շֱ���
r_t_1D = 0.4;    % 1D SBL ʱ�ӷֱ���

% 2D SBL ����
r_v_2D = 0.4;    % 2D SBL �����շֱ���
r_t_2D = 0.4;    % 2D SBL ʱ�ӷֱ���

% �ֲ�1D SBL ����
r_v_coarse_1D = 0.5;    % �ֲ�1D SBL ����������շֱ���
r_t_coarse_1D = 0.5;    % �ֲ�1D SBL ������ʱ�ӷֱ���
r_v_fine_1D = 0.2;      % �ֲ�1D SBL ϸ��������շֱ���
r_t_fine_1D = 0.2;      % �ֲ�1D SBL ϸ����ʱ�ӷֱ���
threshold_ratio_1D = 0.1; % �ֲ�1D SBL ϵ��ѡ����ֵ

% �ֲ�2D SBL ����
r_v_coarse_2D = 0.5;    % �ֲ�2D SBL ����������շֱ���
r_t_coarse_2D = 0.5;    % �ֲ�2D SBL ������ʱ�ӷֱ���
r_v_fine_2D = 0.2;      % �ֲ�2D SBL ϸ��������շֱ���
r_t_fine_2D = 0.2;      % �ֲ�2D SBL ϸ����ʱ�ӷֱ���
threshold_ratio_2D = 0.1; % �ֲ�2D SBL ϵ��ѡ����ֵ

% OMP ����
r_v_OMP = 0.4;    % OMP �����շֱ���
r_t_OMP = 0.4;    % OMP ʱ�ӷֱ���

%% �̶���������
% �̶�SNR
SNR_dB_fixed = 10;  % �̶�SNRΪ15dB
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

%% �������
N_fram = 20;  % ����֡��

%% ��ʼ������洢
NMSE_count_1D = zeros(num_P_values, 1);
NMSE_count_2D = zeros(num_P_values, 1);
NMSE_count_hierarchical_1D = zeros(num_P_values, 1);
NMSE_count_hierarchical_2D = zeros(num_P_values, 1);
NMSE_count_OMP = zeros(num_P_values, 1);
NMSE_count_traditional = zeros(num_P_values, 1);
NMSE_CRLB = zeros(num_P_values, 1);

% ����ʱ��洢
time_1D_SBL = zeros(num_P_values, 1);
time_2D_SBL = zeros(num_P_values, 1);
time_hierarchical_1D = zeros(num_P_values, 1);
time_hierarchical_2D = zeros(num_P_values, 1);
time_OMP = zeros(num_P_values, 1);
time_traditional = zeros(num_P_values, 1);

fprintf('\n=== �㷨������·������P�仯���� ===\n');
fprintf('�̶�SNR: %d dB\n', SNR_dB_fixed);
fprintf('���������N=%d, M=%d, N_fram=%d\n', N, M, N_fram);
fprintf('����·��������Χ: %d ~ %d\n', min(P_values), max(P_values));
fprintf('\n');

%% ��ѭ����������ͬ·������
for P_idx = 1:num_P_values
    P_current = P_values(P_idx);
    
    fprintf('\n--- ����·������: %d ---\n', P_current);
    
    % ��ʼ��ʱ���ۼ���
    total_time_1D = 0;
    total_time_2D = 0;
    total_time_hierarchical_1D = 0;
    total_time_hierarchical_2D = 0;
    total_time_OMP = 0;
    total_time_traditional = 0;
    
    %% Monte Carlo����
    for ifram = 1:N_fram
        if mod(ifram, 5) == 0
            fprintf('  Frame %d/%d\n', ifram, N_fram);
        end
        
        %% ����ŵ���ʼ����ʹ�õ�ǰ·��������
        v_c_init = unifrnd(-v_max, v_max, P_current, 1);
        t_c_init = unifrnd(0, t_max, P_current, 1);
        l_ti = t_c_init .* (M * delta_f);
        q_l_t = exp(-0.1 .* l_ti) ./ sum(exp(-0.1 .* l_ti));
        h_c_init = normrnd(0, q_l_t);
        k_v_init = v_c_init .* (N * T);
        l_t_init = t_c_init .* (M * delta_f);
        
        pow_prof = (1/P_current) * (ones(1, P_current));
        chan_coef = zeros(1, P_current);
        
        %% ��������������
        trans_info_bit = randi([0,1], N_syms_perfram*M_bits, 1);
        data = qammod(reshape(trans_info_bit, M_bits, N_syms_perfram), M_mod, 'gray', 'InputType', 'bit');
        data = data ./ (eng_sqrt);
        x = Generate_2D_data_grid_CE(N, M, data, data_grid);
        x(x_kp-floor(Kp/2):x_kp+floor(Kp/2), x_lp-floor(Lp/2):x_lp+floor(Lp/2)) = sqrt(1000)/(Kp*Lp);
        
        % �����������
        phi_sys = zeros(M*N, P_current);
        kk = 1:N;
        ll = 1:M;
        
        for pp = 1:P_current
            for k = 0:(N-1)
                for l = 1:M
                    v = Sampling_Function_v(N, k, kk-1, k_v_init(pp));
                    t = Sampling_Function_t(M, l-1, ll-1, l_t_init(pp));
                    phi_sys(k*M+l, pp) = sum(sum(v.' * t .* x));
                end
            end
        end
        
        %% �ŵ����
        noise_gen_re = sigma_2_fixed * randn(M*N, 1);
        noise_gen_im = sigma_2_fixed * randn(M*N, 1);
        noise_gen = noise_gen_re + 1i .* noise_gen_im;
        h_c_init_hat = h_c_init .* exp(-2i * pi/M/N * (k_v_init .* l_t_init));
        r = phi_sys * h_c_init + noise_gen;
        y = reshape(r, M, N).';
        
        %% ��ȡ�ضϲ�������
        y_trunc = y(floor(N/2)-floor(Kp/2)-k_max:floor(N/2)+floor(Kp/2)+k_max, floor(M/2)-floor(Lp/2):floor(M/2)+floor(Lp/2)+l_max);
        N_T = 2*k_max + Kp;
        M_T = l_max + Lp;
        y_T = reshape(y_trunc.', N_T*M_T, 1);
        
        %% �㷨1: 1D SBL Channel Estimation
        tic;
        try
            [h_hat_1D, k_v_hat_1D, l_t_hat_1D, virtual_size_1D, Phi_1D, delta_a] = ...
                CE_1D_SBL(sqrt(1000)/(Kp*Lp), k_max+Kp, Lp, M, N, N_T, M_T, y_T, r_v_1D, r_t_1D, k_max, l_max, on_flag);
        catch ME
            fprintf('    1D SBL�㷨����: %s\n', ME.message);
            h_hat_1D = zeros(P_current, 1);
            k_v_hat_1D = zeros(P_current, 1);
            l_t_hat_1D = zeros(P_current, 1);
        end
        total_time_1D = total_time_1D + toc;
        
        %% �㷨2: 2D SBL Channel Estimation
        tic;
        try
            [H_opt_2D, k_v_opt_2D, l_t_opt_2D, virtual_size_2D] = ...
                CE_2D_SBL(sqrt(1000)/(2*Kp*Lp), k_max+Kp, 2, Lp, M, N, N_T, M_T, y_trunc, r_v_2D, r_t_2D, k_max, l_max, on_flag);
        catch ME
            fprintf('    2D SBL�㷨����: %s\n', ME.message);
            H_opt_2D = zeros(P_current, 1);
            k_v_opt_2D = zeros(P_current, 1);
            l_t_opt_2D = zeros(P_current, 1);
        end
        total_time_2D = total_time_2D + toc;
        
        %% �㷨3: �ֲ�1D SBL Channel Estimation
        tic;
        try
            % ���ɴ�����
            k_v_coarse = -k_max:r_v_coarse_1D:k_max;
            l_t_coarse = 0:r_t_coarse_1D:l_max;
            [K_V_coarse, L_T_coarse] = meshgrid(k_v_coarse, l_t_coarse);
            k_v_grid_coarse = K_V_coarse(:);
            l_t_grid_coarse = L_T_coarse(:);
            
            [h_hat_hierarchical_1D, k_v_hat_hierarchical_1D, l_t_hat_hierarchical_1D] = ...
                hierarchical_SBL_refined_1D(sqrt(1000)/(Kp*Lp), k_max+Kp, Lp, M, N, N_T, M_T, y_T, ...
                                           k_v_grid_coarse, l_t_grid_coarse, on_flag);
        catch ME
            fprintf('    �ֲ�1D SBL�㷨����: %s\n', ME.message);
            h_hat_hierarchical_1D = zeros(P_current, 1);
            k_v_hat_hierarchical_1D = zeros(P_current, 1);
            l_t_hat_hierarchical_1D = zeros(P_current, 1);
        end
        total_time_hierarchical_1D = total_time_hierarchical_1D + toc;
        
        %% �㷨4: �ֲ�2D SBL Channel Estimation
        tic;
        try
            % ���ɴ�����
            k_v_coarse_2d = -k_max:r_v_coarse_2D:k_max;
            l_t_coarse_2d = 0:r_t_coarse_2D:l_max;
            [K_V_coarse_2d, L_T_coarse_2d] = meshgrid(k_v_coarse_2d, l_t_coarse_2d);
            k_v_grid_coarse_2d = K_V_coarse_2d(:);
            l_t_grid_coarse_2d = L_T_coarse_2d(:);
            
            [H_opt_hierarchical_2D, k_v_opt_hierarchical_2D, l_t_opt_hierarchical_2D] = ...
                hierarchical_SBL_refined_2D(sqrt(1000)/(2*Kp*Lp), k_max+Kp, 2, Lp, M, N, N_T, M_T, y_trunc, ...
                                           k_v_grid_coarse_2d, l_t_grid_coarse_2d, r_v_coarse_2D, r_t_coarse_2D, ...
                                           k_max, l_max, on_flag);
        catch ME
            fprintf('    �ֲ�2D SBL�㷨����: %s\n', ME.message);
            H_opt_hierarchical_2D = zeros(P_current, 1);
            k_v_opt_hierarchical_2D = zeros(P_current, 1);
            l_t_opt_hierarchical_2D = zeros(P_current, 1);
        end
        total_time_hierarchical_2D = total_time_hierarchical_2D + toc;
        
        %% �㷨5: OMP Channel Estimation
        tic;
        try
            [h_hat_omp, omp_index, k_v_omp, l_t_omp] = ...
                OMP(sqrt(1000)/(Kp*Lp), k_max+Kp, Lp, M, N, N_T, M_T, y_T, r_v_OMP, r_t_OMP, k_max, l_max);
        catch ME
            fprintf('    OMP�㷨����: %s\n', ME.message);
            h_hat_omp = zeros(P_current, 1);
            k_v_omp = zeros(P_current, 1);
            l_t_omp = zeros(P_current, 1);
        end
        total_time_OMP = total_time_OMP + toc;
        
        %% �㷨6: Traditional Impulse Response Estimation
        tic;
        try
            [h_trd, k_v_trd, l_t_trd] = ...
                traditional_impulse(sqrt(1000)/(Kp*Lp), y_trunc, k_max, l_max, sigma_2_p_fixed);
        catch ME
            fprintf('    Traditional�㷨����: %s\n', ME.message);
            h_trd = zeros(P_current, 1);
            k_v_trd = zeros(P_current, 1);
            l_t_trd = zeros(P_current, 1);
        end
        total_time_traditional = total_time_traditional + toc;
        
        %% ����NMSE
        NMSE_nume_1D = 0;
        NMSE_nume_2D = 0;
        NMSE_nume_hierarchical_1D = 0;
        NMSE_nume_hierarchical_2D = 0;
        NMSE_nume_OMP = 0;
        NMSE_nume_traditional = 0;
        NMSE_deno = 0;
        
        for kk = 0:(N-1)
            for ll = 1:M
                % ��ʵ�ŵ���Ӧ
                h_w = sum(Sampling_Function_v(N, kk+1, 1, k_v_init) .* ...
                         Sampling_Function_t(M, ll, 1, l_t_init) .* h_c_init .* ...
                         exp(-2i*pi .* k_v_init .* l_t_init / N / M));
                
                % ���㷨�Ĺ���ֵ
                h_w_hat_1D = sum(Sampling_Function_v(N, kk+1, 1, k_v_hat_1D) .* ...
                                Sampling_Function_t(M, ll, 1, l_t_hat_1D) .* h_hat_1D .* ...
                                exp(-2i*pi .* k_v_hat_1D .* l_t_hat_1D / N / M));
                
                h_w_hat_2D = sum(Sampling_Function_v(N, kk+1, 1, k_v_opt_2D) .* ...
                                Sampling_Function_t(M, ll, 1, l_t_opt_2D) .* H_opt_2D .* ...
                                exp(-2i*pi .* k_v_opt_2D .* l_t_opt_2D / N / M));
                
                h_w_hat_hierarchical_1D = sum(Sampling_Function_v(N, kk+1, 1, k_v_hat_hierarchical_1D) .* ...
                                              Sampling_Function_t(M, ll, 1, l_t_hat_hierarchical_1D) .* h_hat_hierarchical_1D .* ...
                                              exp(-2i*pi .* k_v_hat_hierarchical_1D .* l_t_hat_hierarchical_1D / N / M));
                
                h_w_hat_hierarchical_2D = sum(Sampling_Function_v(N, kk+1, 1, k_v_opt_hierarchical_2D) .* ...
                                              Sampling_Function_t(M, ll, 1, l_t_opt_hierarchical_2D) .* H_opt_hierarchical_2D .* ...
                                              exp(-2i*pi .* k_v_opt_hierarchical_2D .* l_t_opt_hierarchical_2D / N / M));
                
                h_w_omp = sum(Sampling_Function_v(N, kk+1, 1, k_v_omp) .* ...
                             Sampling_Function_t(M, ll, 1, l_t_omp) .* h_hat_omp .* ...
                             exp(-2i*pi .* k_v_omp .* l_t_omp / N / M));
                
                h_w_hat_trd = sum(Sampling_Function_v(N, kk+1, 1, k_v_trd) .* ...
                                 Sampling_Function_t(M, ll, 1, l_t_trd) .* h_trd .* ...
                                 exp(-2i*pi .* k_v_trd .* l_t_trd / N / M));
                
                NMSE_nume_1D = NMSE_nume_1D + abs(h_w - h_w_hat_1D).^2;
                NMSE_nume_2D = NMSE_nume_2D + abs(h_w - h_w_hat_2D).^2;
                NMSE_nume_hierarchical_1D = NMSE_nume_hierarchical_1D + abs(h_w - h_w_hat_hierarchical_1D).^2;
                NMSE_nume_hierarchical_2D = NMSE_nume_hierarchical_2D + abs(h_w - h_w_hat_hierarchical_2D).^2;
                NMSE_nume_OMP = NMSE_nume_OMP + abs(h_w - h_w_omp).^2;
                NMSE_nume_traditional = NMSE_nume_traditional + abs(h_w - h_w_hat_trd).^2;
                NMSE_deno = NMSE_deno + abs(h_w)^2;
            end
        end
        
        NMSE_count_1D(P_idx) = NMSE_count_1D(P_idx) + NMSE_nume_1D / (NMSE_deno * N_fram);
        NMSE_count_2D(P_idx) = NMSE_count_2D(P_idx) + NMSE_nume_2D / (NMSE_deno * N_fram);
        NMSE_count_hierarchical_1D(P_idx) = NMSE_count_hierarchical_1D(P_idx) + NMSE_nume_hierarchical_1D / (NMSE_deno * N_fram);
        NMSE_count_hierarchical_2D(P_idx) = NMSE_count_hierarchical_2D(P_idx) + NMSE_nume_hierarchical_2D / (NMSE_deno * N_fram);
        NMSE_count_OMP(P_idx) = NMSE_count_OMP(P_idx) + NMSE_nume_OMP / (NMSE_deno * N_fram);
        NMSE_count_traditional(P_idx) = NMSE_count_traditional(P_idx) + NMSE_nume_traditional / (NMSE_deno * N_fram);
    end
    
    % ��¼ƽ������ʱ��
    time_1D_SBL(P_idx) = total_time_1D / N_fram;
    time_2D_SBL(P_idx) = total_time_2D / N_fram;
    time_hierarchical_1D(P_idx) = total_time_hierarchical_1D / N_fram;
    time_hierarchical_2D(P_idx) = total_time_hierarchical_2D / N_fram;
    time_OMP(P_idx) = total_time_OMP / N_fram;
    time_traditional(P_idx) = total_time_traditional / N_fram;
    
    fprintf('·������ %d ��ɣ�ƽ������ʱ�䣺1D SBL=%.4fs, 2D SBL=%.4fs, �ֲ�1D=%.4fs, �ֲ�2D=%.4fs, OMP=%.4fs, Traditional=%.4fs\n', ...
        P_current, time_1D_SBL(P_idx), time_2D_SBL(P_idx), time_hierarchical_1D(P_idx), ...
        time_hierarchical_2D(P_idx), time_OMP(P_idx), time_traditional(P_idx));
end

%% �������ͻ�ͼ
NMSE_dB_1D = 10 * log10(NMSE_count_1D);
NMSE_dB_2D = 10 * log10(NMSE_count_2D);
NMSE_dB_hierarchical_1D = 10 * log10(NMSE_count_hierarchical_1D);
NMSE_dB_hierarchical_2D = 10 * log10(NMSE_count_hierarchical_2D);
NMSE_dB_OMP = 10 * log10(NMSE_count_OMP);
NMSE_dB_traditional = 10 * log10(NMSE_count_traditional);

% ��ɫ����
colors = [0,107,182;     % ��ɫ - 1D SBL
          118,41,133;    % ��ɫ - 2D SBL
          234,174,31;    % ��ɫ - �ֲ�1D SBL
          215,94,59;     % ��ɫ - �ֲ�2D SBL
          184,125,182;   % ����ɫ - OMP
          71,90,40;      % ��ɫ - Traditional
          161,27,30]/256; % ��ɫ - ����

%% ����NMSE���ܶԱ�ͼ
figure('Position', [100, 100, 1200, 800]);

plot(P_values, NMSE_dB_1D, '-s', 'Color', colors(1,:), 'LineWidth', 2, 'MarkerSize', 8);
hold on;
plot(P_values, NMSE_dB_2D, '-o', 'Color', colors(2,:), 'LineWidth', 2, 'MarkerSize', 8);
plot(P_values, NMSE_dB_hierarchical_1D, '-^', 'Color', colors(3,:), 'LineWidth', 2, 'MarkerSize', 8);
plot(P_values, NMSE_dB_hierarchical_2D, '-v', 'Color', colors(4,:), 'LineWidth', 2, 'MarkerSize', 8);
plot(P_values, NMSE_dB_OMP, '-d', 'Color', colors(5,:), 'LineWidth', 2, 'MarkerSize', 8);
plot(P_values, NMSE_dB_traditional, '-p', 'Color', colors(6,:), 'LineWidth', 2, 'MarkerSize', 8);
grid on;
xlabel('·������ P', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('NMSE (dB)', 'FontSize', 12, 'FontWeight', 'bold');
title('�㷨NMSE������·�������仯', 'FontSize', 14);
legend('1D SBL', '2D SBL', '�ֲ�1D SBL', '�ֲ�2D SBL', 'OMP', 'Traditional', 'Location', 'best');
set(gca, 'FontSize', 11);