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
r_v = 0.5;
r_t = 0.5;
P = 4;
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

N_fram = 3;
NMSE_count_2D = zeros(length(SNR_dB),1);
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
        
        %% 2D SBL Channel Estimation (Only method used)
        [H_opt_2D,k_v_opt_2D,l_t_opt_2D,virtual_size_2D] = CE_2D_SBL(sqrt(1000)/(2*Kp*Lp),k_max+Kp,2,Lp,M,N,N_T,M_T,y_trunc,r_v,r_t,k_max,l_max,on_flag);

       
        %% count NMSE for 2D SBL only
        NMSE_nume_2D = 0;  
        NMSE_deno = 0;
        for kk = 0:(N-1)
            for ll = 1:M
                if (test==1)
                    h_w = sum(Sampling_Function_v(N,kk+1,1,k_v_test) .* Sampling_Function_t(M,ll,1,l_t_test) .* h_test .* exp(-2i*pi.*k_v_test.*l_t_test/N/M));
                else
                    h_w = sum(Sampling_Function_v(N,kk+1,1,k_v_init) .* Sampling_Function_t(M,ll,1,l_t_init) .* h_c_init .* exp(-2i*pi.*k_v_init.*l_t_init/N/M));
                end
                h_w_hat_2D = sum(Sampling_Function_v(N,kk+1,1,k_v_opt_2D) .* Sampling_Function_t(M,ll,1,l_t_opt_2D) .* H_opt_2D .* exp(-2i*pi.*k_v_opt_2D.*l_t_opt_2D/N/M));
                NMSE_nume_2D = NMSE_nume_2D + abs(h_w - h_w_hat_2D).^2;
                NMSE_deno = NMSE_deno + abs(h_w)^2;
            end 
        end
        
        % No CRLB calculation needed
        
        NMSE_count_2D(iesn0) = NMSE_count_2D(iesn0) + NMSE_nume_2D / (NMSE_deno * N_fram);
        
        display(ifram,'ifram');
        display(iesn0, 'iesn0');
    end
end

%% Results plotting
NMSE_dB_2D = 10 * log10(NMSE_count_2D);

colors = [0,107,182; %?1
          118,41,133; %?2
          234,174,31; %?3
          215,94,59; %?4
          184,125,182; %??5
          71,90,40; %?6
          161,27,30]; %?7
colors = colors/256;

figure()
plot(SNR_dB,NMSE_dB_2D,'-*', 'Color', colors(2,:),'LineWidth',1.5,'MarkerSize',8);
grid on

legend('2D off-grid SBL');
xlabel('SNR(dB)')
ylabel('NMSE(dB)')
title('2D Sparse Bayesian Learning Channel Estimation Performance')

fprintf('\n=== 2D SBL Channel Estimation Results ===\n');
for i = 1:length(SNR_dB)
    fprintf('SNR = %d dB: NMSE = %.2f dB\n', SNR_dB(i), NMSE_dB_2D(i));
end

%% Additional analysis for 2D SBL
fprintf('\n=== 2D SBL Algorithm Details ===\n');
fprintf('Virtual grid size: %d x %d\n', ceil(2 * k_max / r_v), ceil(l_max / r_t));
fprintf('Total virtual points: %d\n', ceil(2 * k_max / r_v) * ceil(l_max / r_t));
fprintf('Doppler resolution: %.2f\n', r_v);
fprintf('Delay resolution: %.2f\n', r_t);
fprintf('Number of paths: %d\n', P);
fprintf('On-grid flag: %d\n', on_flag);

% Results display only, no saving