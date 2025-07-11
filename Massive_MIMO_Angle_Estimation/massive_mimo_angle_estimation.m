close all
clear all
clc

%% Massive MIMO ����ǹ����㷨�ȽϽű�
% ģ��ODDM�ŵ�������ʱ�ӺͶ�����Ƶ�Ƶ�����
% �÷��䷽���(AoD)�ͽ��շ����(AoA)���ʱ�ӺͶ�����Ƶ��
% ������1D SBL, 2D SBL, �ֲ�1D SBL, �ֲ�2D SBL, OMP, Traditional����

%% �㷨�ֱ��ʲ������� (�ɶ�������)
% 1D SBL ����
r_aod_1D = 0.8;    % 1D SBL ����Ƿֱ��� (����)
r_aoa_1D = 0.8;    % 1D SBL ���սǷֱ��� (����)

% 2D SBL ����
r_aod_2D = 0.8;    % 2D SBL ����Ƿֱ���
r_aoa_2D = 0.8;    % 2D SBL ���սǷֱ���

% �ֲ�1D SBL ����
r_aod_coarse_1D = 0.8;    % �ֲ�1D SBL ��������Ƿֱ���
r_aoa_coarse_1D = 0.8;    % �ֲ�1D SBL ��������սǷֱ���
r_aod_fine_1D = 0.4;      % �ֲ�1D SBL ϸ������Ƿֱ���
r_aoa_fine_1D = 0.4;      % �ֲ�1D SBL ϸ������սǷֱ���
threshold_ratio_1D = 0.1;  % �ֲ�1D SBL ϵ��ѡ����ֵ

% �ֲ�2D SBL ����
r_aod_coarse_2D = 0.8;    % �ֲ�2D SBL ��������Ƿֱ���
r_aoa_coarse_2D = 0.8;    % �ֲ�2D SBL ��������սǷֱ���
r_aod_fine_2D = 0.4;      % �ֲ�2D SBL ϸ������Ƿֱ���
r_aoa_fine_2D = 0.4;      % �ֲ�2D SBL ϸ������սǷֱ���
threshold_ratio_2D = 0.1;  % �ֲ�2D SBL ϵ��ѡ����ֵ

% OMP ����
r_aod_OMP = 0.8;    % OMP ����Ƿֱ���
r_aoa_OMP = 0.8;    % OMP ���սǷֱ���

%% Test Mode %%%%%%
test = 0;    % set to 1 when testing

%% Massive MIMO ϵͳ����%%%%%%%%%%
% Nt: ����������
Nt = 32;
% Nr: ����������
Nr = 32;
% ���߼�� (�Բ���Ϊ��λ)
d_lambda = 0.5;
% �ز�Ƶ��
car_freq = 5e9;  % 28 GHz
c = 3e8;          % ����
lambda = c / car_freq;
d = d_lambda * lambda;

% �Ƕȷ�Χ����
aod_max = pi/3;   % ��������ֵ (60��)
aoa_max = pi/3;   % ���ս����ֵ (60��)
P = 5;            % ·����
on_flag = 0;      % set to 0 when on-grid

%% �����ò���
aod_test = [0, 0, 0, 0, 0]';  % �����÷����
aoa_test = [0, 0, 0, 0, 0]';  % �����ý��ս�
h_test = [1, 0, 0, 0, 0]';    % �������ŵ�ϵ��

%% ��Ƶ��������
Kp = 1;
Lp = 1;
x_kp = floor(Nt/2);
x_lp = floor(Nr/2);
% ����λ������
data_grid = ones(Nt, Nr);
data_grid(x_kp-floor(Kp/2)-2:x_kp+floor(Kp/2)+2, x_lp-floor(Lp/2)-2:x_lp+floor(Lp/2)+2) = 0;
% ÿ֡������
N_syms_perfram = sum(sum(data_grid));

% SNR ����������
SNR_dB = 0:5:20;
SNR = 10.^(SNR_dB/10);
sigma_2 = 0.5 ./ SNR;
SNR_p = SNR * 1e5;
sigma_2_p = sigma_2 / 1e5;

%% ��ʼ����������������
N_fram = 10;  % ����֡��
NMSE_count_1D = zeros(length(SNR_dB), 1);
NMSE_count_2D = zeros(length(SNR_dB), 1);
NMSE_count_hierarchical_1D = zeros(length(SNR_dB), 1);
NMSE_count_hierarchical_2D = zeros(length(SNR_dB), 1);
NMSE_count_OMP = zeros(length(SNR_dB), 1);
NMSE_count_traditional = zeros(length(SNR_dB), 1);
NMSE_CRLB = zeros(length(SNR_dB), 1);  % CRLB�����½�

fprintf('\n=== Massive MIMO ����ǹ����㷨�ȽϿ�ʼ ===\n');
fprintf('���������Nt=%d, Nr=%d, P=%d, N_fram=%d\n', Nt, Nr, P, N_fram);
fprintf('\n�㷨�ֱ��ʲ�����\n');
fprintf('1D SBL: r_aod=%.3f, r_aoa=%.3f\n', r_aod_1D, r_aoa_1D);
fprintf('2D SBL: r_aod=%.3f, r_aoa=%.3f\n', r_aod_2D, r_aoa_2D);
fprintf('�ֲ�1D SBL: ������(r_aod=%.3f, r_aoa=%.3f), ϸ����(r_aod=%.3f, r_aoa=%.3f)\n', ...
    r_aod_coarse_1D, r_aoa_coarse_1D, r_aod_fine_1D, r_aoa_fine_1D);
fprintf('�ֲ�2D SBL: ������(r_aod=%.3f, r_aoa=%.3f), ϸ����(r_aod=%.3f, r_aoa=%.3f)\n', ...
    r_aod_coarse_2D, r_aoa_coarse_2D, r_aod_fine_2D, r_aoa_fine_2D);
fprintf('OMP: r_aod=%.3f, r_aoa=%.3f\n', r_aod_OMP, r_aoa_OMP);
fprintf('Traditional: ��ͳ�������ι���\n');
fprintf('\n');

for ifram = 1:N_fram
    
    %% ����ŵ���ʼ��
    aod_c_init = unifrnd(-aod_max, aod_max, P, 1);   % �����
    aoa_c_init = unifrnd(-aoa_max, aoa_max, P, 1);   % ���ս�
    
    % ·������ (ģ��ODDM�еĹ��ʷֲ�)
    path_gains = exp(-0.1 * (1:P)');
    path_gains = path_gains / sum(path_gains);
    h_c_init = normrnd(0, path_gains);  % ����˹�ŵ�ϵ��
    
    for iesn0 = 1:length(SNR_dB)
        fprintf('Frame %d/%d, SNR %d dB\n', ifram, N_fram, SNR_dB(iesn0));
        
        %% ���ɵ�Ƶ�ź�
        pilot_power = sqrt(1000) / (Kp * Lp);
        
        % ����������� (ģ��ODDM��phi_sys)
        phi_sys = zeros(Nr * Nt, P);
        
        for pp = 1:P
            % ����������Ӧ����
            at = exp(1j * 2 * pi * d / lambda * (0:Nt-1)' * sin(aod_c_init(pp)));
            % ����������Ӧ����
            ar = exp(1j * 2 * pi * d / lambda * (0:Nr-1)' * sin(aoa_c_init(pp)));
            
            % ������������ (Kronecker����ʽ)
            phi_sys(:, pp) = kron(at, ar);
        end
        
        %% �ŵ����
        noise_gen_re = sigma_2(iesn0) * randn(Nr * Nt, 1);
        noise_gen_im = sigma_2(iesn0) * randn(Nr * Nt, 1);
        noise_gen = noise_gen_re + 1j * noise_gen_im;
        
        if (test == 1)
            r = phi_sys * h_test;
        else
            r = phi_sys * h_c_init + noise_gen;
        end
        
        y = reshape(r, Nr, Nt).';
        
        %% ��ȡ�ضϲ�������
        aod_range = 4;  % �Ƕ�������Χ (��ӦODDM�е�k_max)
        aoa_range = 4;  % �Ƕ�������Χ (��ӦODDM�е�l_max)
        
        y_trunc = y(floor(Nt/2)-floor(Kp/2)-aod_range:floor(Nt/2)+floor(Kp/2)+aod_range, ...
                   floor(Nr/2)-floor(Lp/2):floor(Nr/2)+floor(Lp/2)+aoa_range);
        N_T = 2 * aod_range + Kp;
        M_T = aoa_range + Lp;
        y_T = reshape(y_trunc.', N_T * M_T, 1);
        
        %% �㷨1: 1D SBL �Ƕȹ���
        [h_hat_1D, aod_hat_1D, aoa_hat_1D, virtual_size_1D, Phi_1D, delta_a] = ...
            MIMO_CE_1D_SBL(pilot_power, aod_range+Kp, Lp, Nr, Nt, N_T, M_T, y_T, ...
                          r_aod_1D, r_aoa_1D, aod_range, aoa_range, on_flag);
        
        %% �㷨2: 2D SBL �Ƕȹ���
        [H_opt_2D, aod_opt_2D, aoa_opt_2D, virtual_size_2D] = ...
            MIMO_CE_2D_SBL(pilot_power, aod_range+Kp, 2, Lp, Nr, Nt, N_T, M_T, y_trunc, ...
                          r_aod_2D, r_aoa_2D, aod_range, aoa_range, on_flag);
        
        %% �㷨3: �ֲ�1D SBL �Ƕȹ���
        % Step 1: ������SBL����
        [h_hat_coarse_1D, aod_hat_coarse_1D, aoa_hat_coarse_1D, virtual_size_coarse_1D, Phi_coarse_1D, delta_coarse_1D] = ...
            MIMO_CE_1D_SBL(pilot_power, aod_range+Kp, Lp, Nr, Nt, N_T, M_T, y_T, ...
                          r_aod_coarse_1D, r_aoa_coarse_1D, aod_range, aoa_range, on_flag);
        
        % Step 2: ʶ������ϵ������λ��
        [sorted_coeff_1D, sorted_idx_1D] = sort(abs(h_hat_coarse_1D), 'descend');
        num_significant_1D = max(1, floor(threshold_ratio_1D * length(h_hat_coarse_1D)));
        significant_idx_1D = sorted_idx_1D(1:num_significant_1D);
        
        % ��ȡ����ϵ����Ӧ�ĽǶ�ֵ
        aod_significant_1D = aod_hat_coarse_1D(significant_idx_1D);
        aoa_significant_1D = aoa_hat_coarse_1D(significant_idx_1D);
        
        % Step 3: ������λ����Χ����ϸ������
        refined_aod_1D = [];
        refined_aoa_1D = [];
        
        for i = 1:length(aod_significant_1D)
            % ��������ϵ����Χ�ľֲ�����
            aod_center = aod_significant_1D(i);
            aoa_center = aoa_significant_1D(i);
            
            % �ڴ�������Χ����ϸ����
            aod_range_local = aod_center + (-r_aod_coarse_1D/2:r_aod_fine_1D:r_aod_coarse_1D/2);
            aoa_range_local = aoa_center + (-r_aoa_coarse_1D/2:r_aoa_fine_1D:r_aoa_coarse_1D/2);
            
            % ȷ���߽�����Ч��Χ��
            aod_range_local = aod_range_local(aod_range_local >= -aod_max & aod_range_local <= aod_max);
            aoa_range_local = aoa_range_local(aoa_range_local >= -aoa_max & aoa_range_local <= aoa_max);
            
            % Ϊ�����򴴽�����
            [AOD_mesh, AOA_mesh] = meshgrid(aod_range_local, aoa_range_local);
            region_aod = AOD_mesh(:);
            region_aoa = AOA_mesh(:);
            
            refined_aod_1D = [refined_aod_1D; region_aod];
            refined_aoa_1D = [refined_aoa_1D; region_aoa];
        end
        
        % ȥ���ظ���
        if ~isempty(refined_aod_1D)
            [refined_grid_1D, unique_idx_1D] = unique([refined_aod_1D, refined_aoa_1D], 'rows');
            refined_aod_1D = refined_grid_1D(:, 1);
            refined_aoa_1D = refined_grid_1D(:, 2);
        end
        
        % Step 4: ��ϸ������ִ��ϸ����SBL
        if length(refined_aod_1D) > 0
            [h_hat_hierarchical_1D, aod_hat_hierarchical_1D, aoa_hat_hierarchical_1D] = ...
                MIMO_hierarchical_SBL_refined_1D(pilot_power, aod_range+Kp, Lp, Nr, Nt, N_T, M_T, y_T, ...
                                                refined_aod_1D, refined_aoa_1D, on_flag);
        else
            % ���˵���������
            h_hat_hierarchical_1D = h_hat_coarse_1D;
            aod_hat_hierarchical_1D = aod_hat_coarse_1D;
            aoa_hat_hierarchical_1D = aoa_hat_coarse_1D;
        end
        
        %% �㷨4: �ֲ�2D SBL �Ƕȹ���
        % Step 1: ������2D SBL����
        [H_opt_coarse_2D, aod_opt_coarse_2D, aoa_opt_coarse_2D, virtual_size_coarse_2D] = ...
            MIMO_CE_2D_SBL(pilot_power, aod_range+Kp, 2, Lp, Nr, Nt, N_T, M_T, y_trunc, ...
                          r_aod_coarse_2D, r_aoa_coarse_2D, aod_range, aoa_range, on_flag);
        
        % Step 2: ʶ������ϵ������2Dλ��
        coeff_magnitude_2D = abs(H_opt_coarse_2D);
        max_coeff_2D = max(coeff_magnitude_2D);
        significant_indices_2D = find(coeff_magnitude_2D > threshold_ratio_2D * max_coeff_2D);
        
        % ����������ת��Ϊ2D����λ��
        N_aod_coarse_2D = ceil(2 * aod_max / r_aod_coarse_2D);
        M_aoa_coarse_2D = ceil(2 * aoa_max / r_aoa_coarse_2D);
        [row_indices_2D, col_indices_2D] = ind2sub([M_aoa_coarse_2D, N_aod_coarse_2D], significant_indices_2D);
        
        % Step 3: ������ϵ����Χ����ϸ������
        [aod_refined_2D, aoa_refined_2D] = create_refined_grid_2D_MIMO(aod_opt_coarse_2D, aoa_opt_coarse_2D, ...
            significant_indices_2D, N_aod_coarse_2D, M_aoa_coarse_2D, r_aod_coarse_2D, r_aoa_coarse_2D, ...
            r_aod_fine_2D, r_aoa_fine_2D, aod_max, aoa_max);
        
        % Step 4: ��ϸ�������Ͻ���ϸ��2D SBL����
        [H_opt_hierarchical_2D, aod_opt_hierarchical_2D, aoa_opt_hierarchical_2D] = ...
            MIMO_hierarchical_SBL_refined_2D(pilot_power, aod_range+Kp, 2, Lp, Nr, Nt, N_T, M_T, y_trunc, ...
                                            aod_refined_2D, aoa_refined_2D, r_aod_fine_2D, r_aoa_fine_2D, ...
                                            aod_range, aoa_range, on_flag);
        
        %% �㷨5: OMP �Ƕȹ���
        [h_hat_OMP, omp_index, aod_hat_OMP, aoa_hat_OMP] = ...
            MIMO_OMP(pilot_power, aod_range+Kp, Lp, Nr, Nt, N_T, M_T, y_T, ...
                    r_aod_OMP, r_aoa_OMP, aod_range, aoa_range);
        
        %% �㷨6: Traditional �������νǶȹ���
        [h_hat_traditional, aod_hat_traditional, aoa_hat_traditional] = ...
            MIMO_traditional_beamforming(pilot_power, y_trunc, aod_range, aoa_range, sigma_2_p(iesn0));
        
        %% ����CRLB�����½�
        % �����ſɱȾ���J (ģ��ODDM�е�ʵ��)
        if exist('virtual_size_1D', 'var') && exist('Phi_1D', 'var') && exist('delta_a', 'var')
            J = zeros(M_T * N_T, virtual_size_1D);
            for m = 1:M_T
                for n = 1:N_T
                    J((n-1)*M_T+m, :) = (MIMO_Array_Response_Tx(Nt, n, 1, aod_hat_1D) .* ...
                                         MIMO_Array_Response_Rx(Nr, m, 1, aoa_hat_1D) .* ...
                                         exp(-2j*pi*aod_hat_1D.*aoa_hat_1D/Nt/Nr)).';
                end
            end
            % ����Э�������
            C_h = (sigma_2(iesn0)^(-1) * (Phi_1D' * Phi_1D) + delta_a^(-1))^(-1);
            C_alpha = J * C_h * J.';
            NMSE_CRLB(iesn0) = trace(abs(C_alpha));
        else
            NMSE_CRLB(iesn0) = 1e-10;  % Ĭ��ֵ
        end
        
        %% ���������㷨��NMSE
        NMSE_nume_1D = 0;
        NMSE_nume_2D = 0;
        NMSE_nume_hierarchical_1D = 0;
        NMSE_nume_hierarchical_2D = 0;
        NMSE_nume_OMP = 0;
        NMSE_nume_traditional = 0;
        NMSE_deno = 0;
        
        for kk = 0:(Nt-1)
            for ll = 1:Nr
                if (test == 1)
                    h_w = sum(MIMO_Array_Response_Tx(Nt, kk+1, 1, aod_test) .* ...
                             MIMO_Array_Response_Rx(Nr, ll, 1, aoa_test) .* h_test .* ...
                             exp(-2j*pi*aod_test.*aoa_test/Nt/Nr));
                else
                    h_w = sum(MIMO_Array_Response_Tx(Nt, kk+1, 1, aod_c_init) .* ...
                             MIMO_Array_Response_Rx(Nr, ll, 1, aoa_c_init) .* h_c_init .* ...
                             exp(-2j*pi*aod_c_init.*aoa_c_init/Nt/Nr));
                end
                
                % ���㷨�Ĺ���ֵ
                h_w_hat_1D = sum(MIMO_Array_Response_Tx(Nt, kk+1, 1, aod_hat_1D) .* ...
                                MIMO_Array_Response_Rx(Nr, ll, 1, aoa_hat_1D) .* h_hat_1D .* ...
                                exp(-2j*pi*aod_hat_1D.*aoa_hat_1D/Nt/Nr));
                
                h_w_hat_2D = sum(MIMO_Array_Response_Tx(Nt, kk+1, 1, aod_opt_2D) .* ...
                                MIMO_Array_Response_Rx(Nr, ll, 1, aoa_opt_2D) .* H_opt_2D .* ...
                                exp(-2j*pi*aod_opt_2D.*aoa_opt_2D/Nt/Nr));
                
                h_w_hat_hierarchical_1D = sum(MIMO_Array_Response_Tx(Nt, kk+1, 1, aod_hat_hierarchical_1D) .* ...
                                              MIMO_Array_Response_Rx(Nr, ll, 1, aoa_hat_hierarchical_1D) .* h_hat_hierarchical_1D .* ...
                                              exp(-2j*pi*aod_hat_hierarchical_1D.*aoa_hat_hierarchical_1D/Nt/Nr));
                
                h_w_hat_hierarchical_2D = sum(MIMO_Array_Response_Tx(Nt, kk+1, 1, aod_opt_hierarchical_2D) .* ...
                                              MIMO_Array_Response_Rx(Nr, ll, 1, aoa_opt_hierarchical_2D) .* H_opt_hierarchical_2D .* ...
                                              exp(-2j*pi*aod_opt_hierarchical_2D.*aoa_opt_hierarchical_2D/Nt/Nr));
                
                h_w_hat_OMP = sum(MIMO_Array_Response_Tx(Nt, kk+1, 1, aod_hat_OMP) .* ...
                                 MIMO_Array_Response_Rx(Nr, ll, 1, aoa_hat_OMP) .* h_hat_OMP .* ...
                                 exp(-2j*pi*aod_hat_OMP.*aoa_hat_OMP/Nt/Nr));
                
                h_w_hat_traditional = sum(MIMO_Array_Response_Tx(Nt, kk+1, 1, aod_hat_traditional) .* ...
                                         MIMO_Array_Response_Rx(Nr, ll, 1, aoa_hat_traditional) .* h_hat_traditional .* ...
                                         exp(-2j*pi*aod_hat_traditional.*aoa_hat_traditional/Nt/Nr));
                
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

%% �����ͼ
NMSE_dB_1D = 10 * log10(NMSE_count_1D);
NMSE_dB_2D = 10 * log10(NMSE_count_2D);
NMSE_dB_hierarchical_1D = 10 * log10(NMSE_count_hierarchical_1D);
NMSE_dB_hierarchical_2D = 10 * log10(NMSE_count_hierarchical_2D);
NMSE_dB_OMP = 10 * log10(NMSE_count_OMP);
NMSE_dB_traditional = 10 * log10(NMSE_count_traditional);
NMSE_dB_CRLB = 10 * log10(NMSE_CRLB);

colors = [0,107,182; %��1
          118,41,133; %��2
          234,174,31; %��3
          215,94,59; %��4
          184,125,182; %����5
          71,90,40; %��6
          161,27,30]; %��7
colors = colors/256;

figure('Position', [100, 100, 800, 600]);
plot(SNR_dB, NMSE_dB_1D, '-v', 'Color', colors(3,:), 'LineWidth', 2, 'MarkerSize', 8);
hold on
plot(SNR_dB, NMSE_dB_2D, '-*', 'Color', colors(2,:), 'LineWidth', 2, 'MarkerSize', 8);
plot(SNR_dB, NMSE_dB_hierarchical_1D, '-s', 'Color', colors(1,:), 'LineWidth', 2, 'MarkerSize', 8);
plot(SNR_dB, NMSE_dB_hierarchical_2D, '-o', 'Color', colors(4,:), 'LineWidth', 2, 'MarkerSize', 8);
plot(SNR_dB, NMSE_dB_OMP, '-d', 'Color', colors(6,:), 'LineWidth', 2, 'MarkerSize', 8);
plot(SNR_dB, NMSE_dB_traditional, '-x', 'Color', colors(7,:), 'LineWidth', 2, 'MarkerSize', 8);
plot(SNR_dB, NMSE_dB_CRLB, '-^', 'Color', colors(5,:), 'LineWidth', 2, 'MarkerSize', 8);
grid on

legend('1D SBL', '2D SBL', '�ֲ�1D SBL', '�ֲ�2D SBL', 'OMP', 'Traditional Beamforming', 'CRLB�½�', 'Location', 'best');
xlabel('SNR (dB)', 'FontSize', 12);
ylabel('NMSE (dB)', 'FontSize', 12);
title('Massive MIMO ����ǹ����㷨���ܱȽ�', 'FontSize', 14);
set(gca, 'FontSize', 11);

%% �����������ʾ
fprintf('\n=== Massive MIMO ����ǹ����㷨�ȽϽ�� ===\n');
fprintf('SNR(dB)\t1D SBL\t\t2D SBL\t\t�ֲ�1D SBL\t�ֲ�2D SBL\tOMP\t\t\tTraditional\tCRLB�½�\n');
fprintf('------\t------\t\t------\t\t----------\t----------\t----------\t----------\t--------\n');
for i = 1:length(SNR_dB)
    fprintf('%d\t\t%.2f\t\t%.2f\t\t%.2f\t\t%.2f\t\t%.2f\t\t%.2f\t\t%.2f\n', ...
        SNR_dB(i), NMSE_dB_1D(i), NMSE_dB_2D(i), NMSE_dB_hierarchical_1D(i), ...
        NMSE_dB_hierarchical_2D(i), NMSE_dB_OMP(i), NMSE_dB_traditional(i), NMSE_dB_CRLB(i));
end

%% �㷨���Ӷȷ���
fprintf('\n=== �㷨���Ӷȷ��� ===\n');
fprintf('1D SBL������: %d\n', ceil(2 * aod_max / r_aod_1D) * ceil(2 * aoa_max / r_aoa_1D));
fprintf('2D SBL������: %d\n', ceil(2 * aod_max / r_aod_2D) * ceil(2 * aoa_max / r_aoa_2D));
fprintf('�ֲ�1D SBL�����������: %d\n', ceil(2 * aod_max / r_aod_coarse_1D) * ceil(2 * aoa_max / r_aoa_coarse_1D));
fprintf('�ֲ�1D SBLϸ���������: %d\n', ceil(2 * aod_max / r_aod_fine_1D) * ceil(2 * aoa_max / r_aoa_fine_1D));
fprintf('�ֲ�2D SBL�����������: %d\n', ceil(2 * aod_max / r_aod_coarse_2D) * ceil(2 * aoa_max / r_aoa_coarse_2D));
fprintf('�ֲ�2D SBLϸ���������: %d\n', ceil(2 * aod_max / r_aod_fine_2D) * ceil(2 * aoa_max / r_aoa_fine_2D));
fprintf('OMP������: %d\n', ceil(2 * aod_max / r_aod_OMP) * ceil(2 * aoa_max / r_aoa_OMP));
fprintf('Traditional Beamforming������: %d\n', (2 * aod_range + 1) * (2 * aoa_range + 1));

fprintf('\n=== ������� ===\n');

%% ������������

%% �ֲ�2D SBL����ϸ������ (Massive MIMO�汾)
function [aod_refined, aoa_refined] = create_refined_grid_2D_MIMO(aod_coarse, aoa_coarse, significant_indices, N_aod_coarse, M_aoa_coarse, r_aod_coarse, r_aoa_coarse, r_aod_fine, r_aoa_fine, aod_max, aoa_max)
    % ΪMassive MIMO������������ϵ����2Dϸ������
    
    % ��ȡ������ṹ
    [aod_bar_coarse, aoa_bar_coarse] = MIMO_First_Order_Linear_Approximation_2D(N_aod_coarse, M_aoa_coarse, aod_max, r_aod_coarse, r_aoa_coarse);
    
    % ��2D�������ҵ�Ψһ������λ��
    significant_aod = [];
    significant_aoa = [];
    
    for idx = 1:length(significant_indices)
        sig_idx = significant_indices(idx);
        
        % ����������ת��Ϊ2D�±�
        [row_idx, col_idx] = ind2sub([M_aoa_coarse, N_aod_coarse], sig_idx);
        
        % ��ȡ����������λ��
        aod_center = aod_bar_coarse(col_idx);
        aoa_center = aoa_bar_coarse(row_idx);
        
        significant_aod = [significant_aod; aod_center];
        significant_aoa = [significant_aoa; aoa_center];
    end
    
    % ȥ���ظ���
    significant_positions = unique([significant_aod, significant_aoa], 'rows');
    significant_aod = significant_positions(:, 1);
    significant_aoa = significant_positions(:, 2);
    
    % ����ϸ������ά��
    N_aod_refined = ceil(2 * aod_max / r_aod_fine);
    M_aoa_refined = ceil(2 * aoa_max / r_aoa_fine);
    
    % ʹ��MIMO_First_Order_Linear_Approximation_2D�����ṹ��ϸ������
    [aod_refined_full, aoa_refined_full] = MIMO_First_Order_Linear_Approximation_2D(N_aod_refined, M_aoa_refined, aod_max, r_aod_fine, r_aoa_fine);
    
    % Ϊ����MIMO_CE_2D_SBL�ṹ���ݣ�����������ϸ������
    % �㷨��ͨ��SBL�����Զ��۽�����������
    aod_refined = aod_refined_full;
    aoa_refined = aoa_refined_full;
end