close all
clear all
clc
%rng(128)

%% ����ģʽ %%%%%%
test = 0;    % ����ʱ����Ϊ1

%% OTFS����%%%%%%%%%%
% N: ʱ�������
N = 32;
% M: Ƶ�����ز���
M = 32;
% M_mod: QAM����ͼ��С
M_mod = 4;
M_bits = log2(M_mod);
% ÿ�����ݷ��ŵ�ƽ������
eng_sqrt = (M_mod==2)+(M_mod~=2)*sqrt((M_mod-1)/6*(2^2));

%% ʱ��-������������ŷ���

% ʱ���Ƶ����Դ
car_fre = 28*10^9;  % �ز�Ƶ��
delta_f = 625*10^3; % ���ز����: 15 KHz
T = 1/delta_f;     % OTFS֡��һ��ʱ����ŵĳ���ʱ��
k_max = 4;         % ϵͳ����
l_max = 4;
t_max = l_max / M / delta_f;
v_max = k_max / N / T;

%% �ֲ�SBL����
% ���������
r_v_coarse = 0.5;    % �����մֱַ���
r_t_coarse = 0.5;    % ʱ�Ӵֱַ���
% ϸ�������
r_v_fine = 0.2;      % ������ϸ�ֱ���
r_t_fine = 0.2;      % ʱ��ϸ�ֱ���
% ����ϸ����ѡ����ֵ
threshold_ratio = 0.1; % ѡ��ǰ10%�Ĺ���ϵ��
P = 4;
on_flag = 0;      % ��������ʱ����Ϊ0

%% �����ò���
k_v_test = [0,0,0,0,0]';
l_t_test = [0,0,0,0,0]';
h_test = [1,0,0,0,0]';

%% ��Ƶ���ŷ���
Kp = 1;
Lp = 1;
x_kp = floor(N/2);
x_lp = floor(M/2);
% OTFSʱ��-�����������ݷ�����2D�����е�����λ��
data_grid = ones(N,M);
data_grid(x_kp-floor(Kp/2)-2*k_max:x_kp+floor(Kp/2)+2*k_max,x_lp-floor(Lp/2)-l_max:x_lp+floor(Lp/2)+l_max)=0;
% ÿ֡������
N_syms_perfram = sum(sum(data_grid));

% ÿ֡������
N_bits_perfram = N_syms_perfram*M_bits;

 
% ����Ⱥ���������
% SNR = P/\sigma^2; P: ������ĸ���ƽ������
SNR_dB = 0:5:20;
SNR = 10.^(SNR_dB/10);
sigma_2 = 0.5 ./SNR;
SNR_p = SNR * 1e5;
sigma_2_p = sigma_2 /1e5;

%% ��ʼ����������������

N_fram = 10;
NMSE_count_coarse = zeros(length(SNR_dB),1);
NMSE_count_hierarchical = zeros(length(SNR_dB),1);
err_ber = zeros(length(SNR_dB),1);
H_mtx = zeros(M*N,M*N);

for ifram = 1:N_fram
    %%  ����ŵ���ʼ��
    v_c_init = unifrnd(-v_max,v_max,P,1);   
    t_c_init = unifrnd(0,t_max,P,1);
    l_ti = t_c_init.*(M * delta_f);
    q_l_t = exp(-0.1.*l_ti)./sum(exp(-0.1.*l_ti));
    h_c_init = normrnd(0,q_l_t);  % ʹ��normrnd������mvnrnd
    k_v_init = v_c_init .*(N*T);
    l_t_init = t_c_init .*(M*delta_f);
    
    for iesn0 = 1:length(SNR_dB)
        pow_prof = (1/P) * (ones(1,P));
        chan_coef = zeros(1,P);
        
        %% ��������������%%%%%
        trans_info_bit = randi([0,1],N_syms_perfram*M_bits,1);
        %%2D QAM�������� %%%%%%%%
        data=qammod(reshape(trans_info_bit,M_bits,N_syms_perfram), M_mod,'gray','InputType','bit');
        data = data./(eng_sqrt);
        x = Generate_2D_data_grid_CE(N,M,data,data_grid);
        x(x_kp-floor(Kp/2):x_kp+floor(Kp/2),x_lp-floor(Lp/2):x_lp+floor(Lp/2))=sqrt(1000)/(Kp*Lp);
        
        % �����������
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
        
        %% �ŵ����%%%%%
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
        
        %% ��ȡ�ضϲ�������%%%%%
        y_trunc = y(floor(N/2)-floor(Kp/2)-k_max:floor(N/2)+floor(Kp/2)+k_max,floor(M/2)-floor(Lp/2):floor(M/2)+floor(Lp/2)+l_max);
        N_T = 2*k_max+Kp;
        M_T = l_max+Lp;
        y_T = reshape(y_trunc.',N_T*M_T,1);        
        
        %% ����1: ������SBL����
        fprintf('Frame %d, SNR %d dB: Starting coarse grid SBL...\n', ifram, SNR_dB(iesn0));
        [h_hat_coarse, k_v_hat_coarse, l_t_hat_coarse, virtual_size_coarse, Phi_coarse, delta_coarse] = ...
            CE_1D_SBL(sqrt(1000)/(Kp*Lp), k_max+Kp, Lp, M, N, N_T, M_T, y_T, r_v_coarse, r_t_coarse, k_max, l_max, on_flag);
        
        %% ����2: ʶ����Ҫϵ������λ��
        [sorted_coeff, sorted_idx] = sort(abs(h_hat_coarse), 'descend');
        num_significant = max(1, floor(threshold_ratio * length(h_hat_coarse)));
        significant_idx = sorted_idx(1:num_significant);
        
        % ��ȡ��Ҫϵ����Ӧ��k_v��l_tֵ
        k_v_significant = k_v_hat_coarse(significant_idx);
        l_t_significant = l_t_hat_coarse(significant_idx);
        
        fprintf('  Selected %d significant coefficients for refinement\n', num_significant);
        
        %% ����3: ����Ҫλ����Χ����ϸ������
        % ����ÿ����Ҫλ����Χ��ϸ������߽�
        refined_k_v = [];
        refined_l_t = [];
        refined_regions = [];
        
        for i = 1:length(k_v_significant)
            % ������Ҫϵ����Χ�ľֲ�����
            k_center = k_v_significant(i);
            l_center = l_t_significant(i);
            
            % �ڴ�������Χ����ϸ����
            k_range = k_center + (-r_v_coarse/2:r_v_fine:r_v_coarse/2);
            l_range = l_center + (-r_t_coarse/2:r_t_fine:r_t_coarse/2);
            
            % ȷ���߽�����Ч��Χ��
            k_range = k_range(k_range >= -k_max & k_range <= k_max);
            l_range = l_range(l_range >= 0 & l_range <= l_max);
            
            % Ϊ�����򴴽�����
            [K_mesh, L_mesh] = meshgrid(k_range, l_range);
            region_k = K_mesh(:);
            region_l = L_mesh(:);
            
            refined_k_v = [refined_k_v; region_k];
            refined_l_t = [refined_l_t; region_l];
            refined_regions = [refined_regions; i*ones(length(region_k), 1)];
        end
        
        % ȥ���ظ���
        [refined_grid, unique_idx] = unique([refined_k_v, refined_l_t], 'rows');
        refined_k_v = refined_grid(:, 1);
        refined_l_t = refined_grid(:, 2);
        
        fprintf('  Created refined grid with %d points\n', length(refined_k_v));
        
        %% ����4: ��ϸ��������ִ��ϸ����SBL
        if length(refined_k_v) > 0
            % Ϊϸ�����񴴽��Զ����������
            [h_hat_fine, k_v_hat_fine, l_t_hat_fine] = ...
                hierarchical_SBL_refined_1D(sqrt(1000)/(Kp*Lp), k_max+Kp, Lp, M, N, N_T, M_T, y_T, ...
                                        refined_k_v, refined_l_t, on_flag);
            
            fprintf('  Completed fine grid SBL\n');
        else
            % ���˵���������
            h_hat_fine = h_hat_coarse;
            k_v_hat_fine = k_v_hat_coarse;
            l_t_hat_fine = l_t_hat_coarse;
        end
        
        %% ����5: ���������ͷֲ㷽����NMSE
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
                
                % ���������
                h_w_hat_coarse = sum(Sampling_Function_v(N,kk+1,1,k_v_hat_coarse) .* Sampling_Function_t(M,ll,1,l_t_hat_coarse) .* h_hat_coarse .* exp(-2i*pi.*k_v_hat_coarse.*l_t_hat_coarse/N/M));
                
                % �ֲ����
                h_w_hat_hierarchical = sum(Sampling_Function_v(N,kk+1,1,k_v_hat_fine) .* Sampling_Function_t(M,ll,1,l_t_hat_fine) .* h_hat_fine .* exp(-2i*pi.*k_v_hat_fine.*l_t_hat_fine/N/M));
                
                NMSE_nume_coarse = NMSE_nume_coarse + abs(h_w - h_w_hat_coarse).^2;
                NMSE_nume_hierarchical = NMSE_nume_hierarchical + abs(h_w - h_w_hat_hierarchical).^2;
                NMSE_deno = NMSE_deno + abs(h_w)^2;
            end 
        end
        
        % ����ҪCRLB����
        
        NMSE_count_coarse(iesn0) = NMSE_count_coarse(iesn0) + NMSE_nume_coarse / (NMSE_deno * N_fram);
        NMSE_count_hierarchical(iesn0) = NMSE_count_hierarchical(iesn0) + NMSE_nume_hierarchical / (NMSE_deno * N_fram);
        
        fprintf('  Coarse NMSE: %.4f, Hierarchical NMSE: %.4f\n', NMSE_nume_coarse/NMSE_deno, NMSE_nume_hierarchical/NMSE_deno);
        
        display(ifram,'ifram');
        display(iesn0, 'iesn0');
    end
end

%%�����ͼ
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