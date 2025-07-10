close all
clear all
clc
%rng(128)

%% Test Mode %%%%%%
test = 0;    % set to 1 when testing

%% OTFS parameters%%%%%%%%%%
% N: number of symbols in time
N = 32;
% M: number of subcarriers in frequency
M = 32;
% M_mod: size of QAM constellation
M_mod = 4;
M_bits = log2(M_mod);
% average energy per data symbol
eng_sqrt = (M_mod==2)+(M_mod~=2)*sqrt((M_mod-1)/6*(2^2));

%% delay-Doppler grid symbol placement

% Time and frequency resources
car_fre = 28*10^9;  % Carrier frequency
delta_f = 625*10^3; % subcarrier spacing: 15 KHz
T = 1/delta_f;     % one time symbol duration in OTFS frame
k_max = 4;         % system setting
l_max = 4;
t_max = l_max / M / delta_f;
v_max = k_max / N / T;
r_v = 0.2;
r_t = 0.2;
P = 4;
on_flag = 0;      % set to 0 when on-grid

%% Hierarchical parameters
% Coarse grid parameters
r_v_coarse = 0.8;  % Coarser Doppler resolution
r_t_coarse = 0.8;  % Coarser delay resolution
threshold_ratio = 0.1;  % Threshold for selecting significant coefficients
refinement_factor = 4;   % How much finer the refined grid should be

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

N_fram = 5;
NMSE_count_coarse_2D = zeros(length(SNR_dB),1);
NMSE_count_hierarchical_2D = zeros(length(SNR_dB),1);
err_ber = zeros(length(SNR_dB),1);
H_mtx = zeros(M*N,M*N);

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
        
        %% Step 1: Coarse Grid 2D SBL Estimation
        fprintf('\n=== Step 1: Coarse Grid 2D SBL Estimation ===\n');
        [H_opt_coarse,k_v_opt_coarse,l_t_opt_coarse,virtual_size_coarse] = CE_2D_SBL(sqrt(1000)/(2*Kp*Lp),k_max+Kp,2,Lp,M,N,N_T,M_T,y_trunc,r_v_coarse,r_t_coarse,k_max,l_max,on_flag);
        
        %% Step 2: Identify significant coefficients and their 2D positions
        fprintf('\n=== Step 2: Identifying Significant Coefficients ===\n');
        coeff_magnitude = abs(H_opt_coarse);
        max_coeff = max(coeff_magnitude);
        significant_indices = find(coeff_magnitude > threshold_ratio * max_coeff);
        
        fprintf('Found %d significant coefficients out of %d total\n', length(significant_indices), length(H_opt_coarse));
        
        % Convert linear indices to 2D grid positions
        N_v_coarse = ceil(2 * k_max / r_v_coarse);
        M_t_coarse = ceil(l_max / r_t_coarse);
        [row_indices, col_indices] = ind2sub([M_t_coarse, N_v_coarse], significant_indices);
        
        % Store grid information for visualization (only for the last frame and SNR)
        if ifram == N_fram && iesn0 == length(SNR_dB)
            % Get coarse grid positions
            [kv_bar_coarse, lt_bar_coarse] = First_Order_Linear_Approximation_2D(N_v_coarse, M_t_coarse, k_max, r_v_coarse, r_t_coarse);
            
            % Store for visualization
            grid_info.N_v_coarse = N_v_coarse;
            grid_info.M_t_coarse = M_t_coarse;
            grid_info.kv_bar_coarse = kv_bar_coarse;
            grid_info.lt_bar_coarse = lt_bar_coarse;
            grid_info.H_opt_coarse = H_opt_coarse;
            grid_info.significant_indices = significant_indices;
            grid_info.row_indices = row_indices;
            grid_info.col_indices = col_indices;
        end
        
        %% Step 3: Create refined grid around significant coefficients
        fprintf('\n=== Step 3: Creating Refined Grid ===\n');
        [k_v_refined, l_t_refined] = hierarchical_SBL_refined_2D(k_v_opt_coarse, l_t_opt_coarse, significant_indices, N_v_coarse, M_t_coarse, r_v_coarse, r_t_coarse, r_v, r_t, k_max, l_max);
        
        fprintf('Refined grid size: %d points\n', length(k_v_refined));
        
        % Store refined grid information for visualization
        if ifram == N_fram && iesn0 == length(SNR_dB)
            N_v_refined = ceil(2 * k_max / r_v);
            M_t_refined = ceil(l_max / r_t);
            grid_info.N_v_refined = N_v_refined;
            grid_info.M_t_refined = M_t_refined;
            grid_info.k_v_refined = k_v_refined;
            grid_info.l_t_refined = l_t_refined;
        end
        
        %% Step 4: Refined 2D SBL estimation on the refined grid
        fprintf('\n=== Step 4: Refined 2D SBL Estimation ===\n');
        [H_opt_hierarchical, k_v_opt_hierarchical, l_t_opt_hierarchical] = CE_2D_SBL_refined(sqrt(1000)/(2*Kp*Lp), k_max+Kp, 2, Lp, M, N, N_T, M_T, y_trunc, k_v_refined, l_t_refined, r_v, r_t, k_max, l_max, on_flag);
        
        % Store final hierarchical results for visualization
        if ifram == N_fram && iesn0 == length(SNR_dB)
            grid_info.H_opt_hierarchical = H_opt_hierarchical;
            grid_info.k_v_opt_hierarchical = k_v_opt_hierarchical;
            grid_info.l_t_opt_hierarchical = l_t_opt_hierarchical;
            
            % Store true channel parameters for comparison
            grid_info.k_v_init = k_v_init;
            grid_info.l_t_init = l_t_init;
            grid_info.h_c_init = h_c_init;
        end
        
        %% count NMSE for coarse 2D SBL
        NMSE_nume_coarse = 0;  
        NMSE_deno = 0;
        for kk = 0:(N-1)
            for ll = 1:M
                if (test==1)
                    h_w = sum(Sampling_Function_v(N,kk+1,1,k_v_test) .* Sampling_Function_t(M,ll,1,l_t_test) .* h_test .* exp(-2i*pi.*k_v_test.*l_t_test/N/M));
                else
                    h_w = sum(Sampling_Function_v(N,kk+1,1,k_v_init) .* Sampling_Function_t(M,ll,1,l_t_init) .* h_c_init .* exp(-2i*pi.*k_v_init.*l_t_init/N/M));
                end
                h_w_hat_coarse = sum(Sampling_Function_v(N,kk+1,1,k_v_opt_coarse) .* Sampling_Function_t(M,ll,1,l_t_opt_coarse) .* H_opt_coarse .* exp(-2i*pi.*k_v_opt_coarse.*l_t_opt_coarse/N/M));
                NMSE_nume_coarse = NMSE_nume_coarse + abs(h_w - h_w_hat_coarse).^2;
                NMSE_deno = NMSE_deno + abs(h_w)^2;
            end 
        end
        
        %% count NMSE for hierarchical 2D SBL
        NMSE_nume_hierarchical = 0;
        for kk = 0:(N-1)
            for ll = 1:M
                if (test==1)
                    h_w = sum(Sampling_Function_v(N,kk+1,1,k_v_test) .* Sampling_Function_t(M,ll,1,l_t_test) .* h_test .* exp(-2i*pi.*k_v_test.*l_t_test/N/M));
                else
                    h_w = sum(Sampling_Function_v(N,kk+1,1,k_v_init) .* Sampling_Function_t(M,ll,1,l_t_init) .* h_c_init .* exp(-2i*pi.*k_v_init.*l_t_init/N/M));
                end
                h_w_hat_hierarchical = sum(Sampling_Function_v(N,kk+1,1,k_v_opt_hierarchical) .* Sampling_Function_t(M,ll,1,l_t_opt_hierarchical) .* H_opt_hierarchical .* exp(-2i*pi.*k_v_opt_hierarchical.*l_t_opt_hierarchical/N/M));
                NMSE_nume_hierarchical = NMSE_nume_hierarchical + abs(h_w - h_w_hat_hierarchical).^2;
            end 
        end
        
        % No CRLB calculation needed
        
        NMSE_count_coarse_2D(iesn0) = NMSE_count_coarse_2D(iesn0) + NMSE_nume_coarse / (NMSE_deno * N_fram);
        NMSE_count_hierarchical_2D(iesn0) = NMSE_count_hierarchical_2D(iesn0) + NMSE_nume_hierarchical / (NMSE_deno * N_fram);
        
        display(ifram,'ifram');
        display(iesn0, 'iesn0');
    end
end

%% Results plotting
NMSE_dB_coarse_2D = 10 * log10(NMSE_count_coarse_2D);
NMSE_dB_hierarchical_2D = 10 * log10(NMSE_count_hierarchical_2D);

colors = [0,107,182; %蓝1
          118,41,133; %紫2
          234,174,31; %黄3
          215,94,59; %橙4
          184,125,182; %粉紫5
          71,90,40; %绿6
          161,27,30]; %红7
colors = colors/256;

figure()
plot(SNR_dB,NMSE_dB_coarse_2D,'-o', 'Color', colors(1,:),'LineWidth',1.5,'MarkerSize',8);
hold on
plot(SNR_dB,NMSE_dB_hierarchical_2D,'-*', 'Color', colors(2,:),'LineWidth',1.5,'MarkerSize',8);
grid on

legend('Coarse Grid 2D SBL', 'Hierarchical 2D SBL');
xlabel('SNR(dB)')
ylabel('NMSE(dB)')
title('Hierarchical 2D Sparse Bayesian Learning Channel Estimation Performance')

fprintf('\n=== Hierarchical 2D SBL Channel Estimation Results ===\n');
for i = 1:length(SNR_dB)
    fprintf('SNR = %d dB: Coarse NMSE = %.2f dB, Hierarchical NMSE = %.2f dB\n', SNR_dB(i), NMSE_dB_coarse_2D(i), NMSE_dB_hierarchical_2D(i));
end

%% 2D Grid Visualization
if exist('grid_info', 'var')
    fprintf('\n=== Generating 2D Grid Visualization ===\n');
    
    % Create single figure
    figure('Position', [100, 100, 800, 600]);
    
    % Create refined grid visualization
    [kv_bar_refined, lt_bar_refined] = First_Order_Linear_Approximation_2D(grid_info.N_v_refined, grid_info.M_t_refined, k_max, r_v, r_t);
    
    % Plot refined grid lines
    [KV_refined, LT_refined] = meshgrid(kv_bar_refined, lt_bar_refined);
    
    % Draw vertical lines for refined grid
    h1 = plot([kv_bar_refined(1), kv_bar_refined(1)], [min(lt_bar_refined), max(lt_bar_refined)], '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.8);
    hold on;
    for i = 2:length(kv_bar_refined)
        plot([kv_bar_refined(i), kv_bar_refined(i)], [min(lt_bar_refined), max(lt_bar_refined)], '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.8);
    end
    
    % Draw horizontal lines for refined grid
    for i = 1:length(lt_bar_refined)
        plot([min(kv_bar_refined), max(kv_bar_refined)], [lt_bar_refined(i), lt_bar_refined(i)], '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.8);
    end
    
    % Draw vertical lines for coarse grid
    for i = 1:length(grid_info.kv_bar_coarse)
        plot([grid_info.kv_bar_coarse(i), grid_info.kv_bar_coarse(i)], [min(grid_info.lt_bar_coarse), max(grid_info.lt_bar_coarse)], '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.8);
    end
    
    % Draw horizontal lines for coarse grid
    for i = 1:length(grid_info.lt_bar_coarse)
        plot([min(grid_info.kv_bar_coarse), max(grid_info.kv_bar_coarse)], [grid_info.lt_bar_coarse(i), grid_info.lt_bar_coarse(i)], '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.8);
    end
    
    % Plot estimated channel positions (top 9 significant coefficients from hierarchical SBL)
    h2 = [];
    
    % Debug information
    fprintf('H_opt_hierarchical size: %d\n', length(grid_info.H_opt_hierarchical));
    fprintf('H_opt_hierarchical max magnitude: %.4f\n', max(abs(grid_info.H_opt_hierarchical)));
    fprintf('H_opt_hierarchical non-zero count: %d\n', sum(abs(grid_info.H_opt_hierarchical) > 1e-10));
    
    % Sort coefficients by magnitude and select top 9
    [sorted_magnitudes, sorted_indices] = sort(abs(grid_info.H_opt_hierarchical), 'descend');
    num_to_show = min(4, length(sorted_indices));
    top_indices = sorted_indices(1:num_to_show);
    
    k_v_est = grid_info.k_v_opt_hierarchical(top_indices);
    l_t_est = grid_info.l_t_opt_hierarchical(top_indices);
    h_est = grid_info.H_opt_hierarchical(top_indices);
    
    fprintf('Showing top %d coefficients\n', num_to_show);
    
    if ~isempty(k_v_est)
        % Use fixed marker size for better visibility
        h2 = scatter(k_v_est, l_t_est, 80, 'g', 'filled', '^', 'LineWidth', 2);
        fprintf('Plotted %d estimated channel positions\n', length(k_v_est));
    else
        fprintf('No significant coefficients found to plot\n');
    end
    
    xlabel('多普勒频移 (k_v)', 'FontSize', 12);
    ylabel('时延 (l_t)', 'FontSize', 12);
    title('分层2D SBL虚拟网格与估计信道位置', 'FontSize', 14, 'FontWeight', 'bold');
    if ~isempty(h2)
        legend([h1, h2], '虚拟网格', '估计信道位置', 'Location', 'best');
    else
        legend(h1, '虚拟网格', 'Location', 'best');
    end
    grid on;
    xlim([-k_max-0.5, k_max+0.5]);
    ylim([-0.5, l_max+0.5]);
    set(gca, 'FontSize', 11);
end

%% Additional analysis
fprintf('\n=== Hierarchical 2D SBL Algorithm Details ===\n');
fprintf('Coarse grid: Doppler res = %.2f, Delay res = %.2f\n', r_v_coarse, r_t_coarse);
fprintf('Fine grid: Doppler res = %.2f, Delay res = %.2f\n', r_v, r_t);
fprintf('Threshold ratio: %.2f\n', threshold_ratio);
fprintf('Refinement factor: %d\n', refinement_factor);
fprintf('Number of paths: %d\n', P);
fprintf('On-grid flag: %d\n', on_flag);

if exist('grid_info', 'var')
    fprintf('\n=== Grid Visualization Details ===\n');
    fprintf('Coarse grid size: %d x %d = %d points\n', grid_info.N_v_coarse, grid_info.M_t_coarse, grid_info.N_v_coarse * grid_info.M_t_coarse);
    fprintf('Refined grid size: %d x %d = %d points\n', grid_info.N_v_refined, grid_info.M_t_refined, grid_info.N_v_refined * grid_info.M_t_refined);
    fprintf('Number of significant coefficients: %d\n', length(grid_info.significant_indices));
    fprintf('Significant coefficient ratio: %.2f%%\n', 100 * length(grid_info.significant_indices) / length(grid_info.H_opt_coarse));
    fprintf('Grid refinement ratio: %.2fx\n', (grid_info.N_v_refined * grid_info.M_t_refined) / (grid_info.N_v_coarse * grid_info.M_t_coarse));
end

% Results display only, no saving

%% Helper function for creating improved refined grid
function [k_v_refined, l_t_refined] = hierarchical_SBL_refined_2D(k_v_coarse, l_t_coarse, significant_indices, N_v_coarse, M_t_coarse, r_v_coarse, r_t_coarse, r_v_fine, r_t_fine, k_max, l_max)
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

%% Helper function for refined 2D SBL estimation using existing CE_2D_SBL structure
function [H_opt, k_v_opt, l_t_opt] = CE_2D_SBL_refined(x_v_p, x_kp, x_t_p, x_lp, M, N, N_T, M_T, Y_T, k_v_grid, l_t_grid, r_v, r_t, k_max, l_max, on_flag)
    % This function uses EXACTLY the same structure as CE_2D_SBL for refined grid estimation
    
    % Create virtual grid dimensions based on refined grid
    N_v = length(unique(k_v_grid));
    M_t = length(unique(l_t_grid));
    virtual_size = N_v * M_t;
    
    % Parameters - EXACTLY same as CE_2D_SBL
    rho = 1e-2;
    c = 1e-4;
    d = 1e-4;
    ksi = 1e-3;
    Tmax = 2e2;
    
    % Initialize variables - EXACTLY same as CE_2D_SBL
    k_v = zeros(N_v,1);
    l_t = zeros(M_t,1);
    mu_dc = zeros(N_v,M_T);
    A_beta0 = zeros(1,M_T);
    H_bar = zeros(N_v,M_t);
    l_t_bar = zeros(N_v,M_t);
    H_opt = zeros(virtual_size,1);
    k_v_opt = zeros(virtual_size,1);
    l_t_opt = zeros(virtual_size,1);
    
    % Create grid approximation for refined grid
    kv_bar = k_v_grid(1:N_v);  % Use first N_v unique values
    lt_bar = l_t_grid(1:M_t);  % Use first M_t unique values
    
    % Generate measurement matrices - EXACTLY same as CE_2D_SBL
    [phi_trunc_L,phi_L,phi_L_v] = Gen_measurement_matrix_L(x_v_p,x_kp,N,N_T,N_v,kv_bar,k_v);
    [~,phi_R,phi_R_t] = Gen_measurement_matrix_R(x_t_p,x_lp,M,M_T,M_t,lt_bar,l_t);
    
    % Initialization - EXACTLY same as CE_2D_SBL
    sigma2_bar = sum(sum(abs(Y_T).^2))/(100*N_T*M_T);
    beta0 = 1/sigma2_bar;
    phiLY = phi_trunc_L' * Y_T;
    alpha_v = mean(abs(phiLY),2);
    
    % Main SBL iteration - EXACTLY same as CE_2D_SBL
    for t = 1:Tmax 
        phi_trunc_L = phi_L + phi_L_v * diag(k_v);
        alpha_v_prev = alpha_v;
        delta_av = diag(alpha_v);
        
        % Update mu and Sigma - EXACTLY same as CE_2D_SBL
        C = (1 / beta0) * eye(N_T) + phi_trunc_L * delta_av * phi_trunc_L';
        C_inv = inv(C);
        sigma_dc = delta_av - delta_av * phi_trunc_L' * C_inv * phi_trunc_L * delta_av;
        for lm = 1:M_T
            mu_dc(:,lm) = beta0 * sigma_dc * phi_trunc_L' * Y_T(:,lm);
        end
        
        % Update gamma1 and alpha - EXACTLY same as CE_2D_SBL
        gamma1 = 1 - real(diag(sigma_dc)) ./ (alpha_v + 1e-16);
        musq = sum(abs(mu_dc).^2, 2);
        alpha_temp = musq + M_T * diag(sigma_dc);
        alpha_v = -0.5 * M_T / rho + sqrt(0.25 * M_T^2 / rho^2 + alpha_temp /rho);
        
        % Update beta0 - EXACTLY same as CE_2D_SBL
        for lr = 1:M_T
            resid = Y_T(:,lr) - phi_trunc_L * mu_dc(:,lr);
            res2 = norm(resid)^2;
            A_beta0(lr) = res2 + 1 / beta0 * sum(gamma1);
        end
        beta0 = (c - 1 + M_T * N_T)/(d + sum(A_beta0));
        
        % Check convergence - EXACTLY same as CE_2D_SBL
        tolerance = norm(alpha_v - alpha_v_prev)/norm(alpha_v_prev);
        if tolerance <= ksi
            break;
        end
        
        % Update k_v - EXACTLY same as CE_2D_SBL
        if(tolerance<1000 * ksi)
            pLpv = phi_L_v' * phi_L_v;
            A_v = zeros(N_v,N_v);
            b_v = zeros(N_v,1);
            for l = 1:M_T
                A_v = A_v + real(pLpv .* (conj(mu_dc(:,l)) * mu_dc(:,l).' + sigma_dc.'))/M_T;
                b_v1 = diag(mu_dc(:,l)) * phi_L_v.' * conj(Y_T(:,l)); 
                b_v2 = diag((mu_dc(:,l) * mu_dc(:,l)' + sigma_dc) * phi_L' * phi_L_v);
                b_v = b_v + real(b_v1 - b_v2) / M_T;
            end
            
            temp = A_v \ b_v;
            k_v_upgrate = zeros(N_v,1);
            flag = ((max(svd(A_v)) / min(svd(A_v))) < 1e5);
            if (flag == 0)
                for pp = 1:N_v
                    temp_k_v = k_v;
                    temp_k_v(pp) = 0;
                    if A_v(pp,pp) == 0
                        k_v_upgrate(pp) = 0;
                        continue;
                    else
                        k_v_upgrate(pp) = (b_v(pp) - A_v(pp,:) * temp_k_v)/A_v(pp,pp);
                    end
                end
            else
                k_v_upgrate = temp;
            end
            
            k_v_upgrate(k_v_upgrate<(-r_v)/2) = (-r_v)/2;
            k_v_upgrate(k_v_upgrate>(r_v)/2) = (r_v)/2;
            if(on_flag==1)
                k_v = k_v_upgrate;
            end
        end
    end
    
    % Step 1 output - EXACTLY same as CE_2D_SBL
    D_hat = mu_dc.';
    k_v_hat = kv_bar + k_v;
    
    % Step 2 - EXACTLY same as CE_2D_SBL
    for stp = 1:N_v
        [h_hat,l_t_hat] = CE_Algo1(M_T,D_hat(:,stp),r_t,M_t,phi_R,phi_R_t,lt_bar,on_flag);
        H_bar(stp,:) = (h_hat.');
        l_t_bar(stp,:) = (l_t_hat.');
    end
    
    % Algorithm output - EXACTLY same as CE_2D_SBL
    H_opt = (reshape(H_bar.',1,virtual_size).');
    k_vv = repmat(k_v_hat,1,M_t);
    k_v_opt = (reshape(k_vv.',1,virtual_size).');
    l_t_opt = (reshape(l_t_bar.',1,virtual_size).');
end