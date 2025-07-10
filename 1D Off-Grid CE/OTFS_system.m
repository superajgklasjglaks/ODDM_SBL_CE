close all
clear all
rng(128)

%% Test Mode %%%%%%
test = 0;

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
car_fre = 3*10^9;  % Carrier frequency
delta_f = 15*10^3; % subcarrier spacing: 15 KHz
T = 1/delta_f;     % one time symbol duration in OTFS frame
k_max = 3;         % system setting
l_max = 4;
t_max = l_max / M / delta_f;
v_max = k_max / N / T;
r_v = 0.2;
r_t = 0.2;
P = 5;
% r_v = 1;
% r_t = 1;

%% parameters for test use
% k_v_test = [0.1,0,0,0,0]';
% l_t_test = [0.7,0,0,0,0]';
% h_test = [5,0,0,0,0]';

k_v_test = [0,0,0,0,0]';
l_t_test = [0,0,0,0,0]';
h_test = [1,0,0,0,0]';

% k_v_test = [0,0,0,0,0]';
% k_v_test = [0.1,0.3,0.5,0.7,0.9]';
% l_t_test = [0.9,0.7,0.5,0.3,0.1]';
% h_test = [1,2,3,4,5]';

% k_v_test = [-2.14,-1.31,0.52,1.7,2.89]';
% l_t_test = [0.92,1.72,1.45,2.33,0.11]';
% h_test = [1,2,3,4,5]';

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

% data_grid(x_kp-floor(Kp/2):x_kp+floor(Kp/2),x_lp-floor(Lp/2):x_lp+floor(Lp/2))=1000/(Kp*Lp);

% number of bits per frame
N_bits_perfram = N_syms_perfram*M_bits;

 
% SNR and variance of the noise
% SNR = P/\sigma^2; P: avg. power of albhabet transmitted
% SNR_dB = 20:-5:0;
SNR_dB = 0:5:30;
SNR = 10.^(SNR_dB/10);
sigma_2 = 0.5 ./SNR;


%% Initializing simulation error count variables

N_fram = 10;
NMSE_count = zeros(length(SNR_dB),1);
NMSE_dB = zeros(length(SNR_dB),1);
err_ber = zeros(length(SNR_dB),1);
H_mtx = zeros(M*N,M*N);
% Initialize_error_count_variables;

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
%         chan_coef = sqrt(pow_prof).*(sqrt(1/2) * (randn(1,P)+1i*randn(1,P)));
%         current_frame_number(iesn0)=ifram;
        
        %% random input bits generation%%%%%
        trans_info_bit = randi([0,1],N_syms_perfram*M_bits,1);
        %%2D QAM symbols generation %%%%%%%%
        data=qammod(reshape(trans_info_bit,M_bits,N_syms_perfram), M_mod,'gray','InputType','bit');
        data = data./(eng_sqrt);
        x = Generate_2D_data_grid_CE(N,M,data,data_grid);
%         data=qammod(reshape(trans_info_bit,M_bits,N_syms_perfram), M_mod,'gray','InputType','bit');
        x(x_kp-floor(Kp/2):x_kp+floor(Kp/2),x_lp-floor(Lp/2):x_lp+floor(Lp/2))=sqrt(1000)/(Kp*Lp);
        
        % calculate measument matrix
        phi_sys = zeros(M*N,P);
        
        for pp=1:P
            for k = 0:(N-1)
                for l = 1:M
                    for kk = 1:N
                        for ll = 1:M
                            if (test==1)
                                phi_sys(k*M+l,pp) = phi_sys(k*M+l,pp) + x(kk,ll) * Sampling_Function_v(N,k,kk-1,k_v_test(pp)) * Sampling_Function_t(M,l-1,ll-1,l_t_test(pp));
                            else
                                phi_sys(k*M+l,pp) = phi_sys(k*M+l,pp) + x(kk,ll) * Sampling_Function_v(N,k,kk-1,k_v_init(pp)) * Sampling_Function_t(M,l-1,ll-1,l_t_init(pp));
                            end
                        end
                    end
                end
            end
        end
        
        %% channel output%%%%%
%         noise_gen_re = mvnrnd(zeros(M*N,1),sigma_2(iesn0).*eye(M*N))';
%         noise_gen_im = mvnrnd(zeros(M*N,1),sigma_2(iesn0).*eye(M*N))';
        noise_gen_re = sigma_2(iesn0) * randn(M*N,1);
        noise_gen_im = sigma_2(iesn0) * randn(M*N,1);
        noise_gen = noise_gen_re + 1i.* noise_gen_im;
%         noise_gen = noise_gen./10;
        h_test_hat = h_test .* exp(-2i * pi/M/N * (k_v_test.*l_t_test));
        h_c_init_hat = h_c_init .* exp(-2i * pi/M/N * (k_v_init.*l_t_init));
        if (test==1)
            r = phi_sys * h_test;
%             r = phi_sys * h_test + noise_gen;
        else
%             r = phi_sys * h_c_init;
            r = phi_sys * h_c_init + noise_gen;
        end
        y = reshape(r,M,N).';
%         y_data = y .* data_grid;
        
        %% get truncation measument matrix%%%%%
%         y_trunc = zeros(2*k_max+Kp,l_max+Lp);
        y_trunc = y(floor(N/2)-floor(Kp/2)-k_max:floor(N/2)+floor(Kp/2)+k_max,floor(M/2)-floor(Lp/2):floor(M/2)+floor(Lp/2)+l_max);
        N_T = 2*k_max+Kp;
        M_T = l_max+Lp;
        y_T = reshape(y_trunc.',N_T*M_T,1);
        [h_hat,k_v_hat,l_t_hat,virtual_size] = CE_1D_SBL(sqrt(1000)/(Kp*Lp),k_max+Kp,Lp,M,N,N_T,M_T,y_T,r_v,r_t,k_max,l_max);  % paper miss
        
        %% MPA detection
%         [~,idx] = sort(h_hat, 'descend');
%         idx = idx(1:4*P);
%         h_res = h_hat(idx);
%         k_v_res = k_v_hat(idx);
%         l_t_res = l_t_hat(idx);
%         for pp = 1:4*P
%             for ky = 0:N-1
%                 for ly =1:M
%                     for kx = 0:N-1
%                         for lx = 1:M
%                             if (test==1)
%                                 H_mtx(ky * M + ly,kx * M + lx) = H_mtx(ky * M + ly,kx * M + lx) + h_test(pp) * Sampling_Function_v(N,ky,kx,k_v_test(pp)) * Sampling_Function_t(M,ly,lx,l_t_test(pp))* exp(-2i * pi/M/N * (k_v_test(pp)*l_t_test(pp)));
%                             else
%                                 H_mtx(ky * M + ly,kx * M + lx) = H_mtx(ky * M + ly,kx * M + lx) + h_res(pp) * Sampling_Function_v(N,ky,kx,k_v_res(pp)) * Sampling_Function_t(M,ly,lx,l_t_res(pp))* exp(-2i * pi/M/N * (k_v_res(pp)*l_t_res(pp)));
%                             end
%                         end
%                     end
%                 end
%             end
%         end
%         
%         x_mtx = reshape(x.',M*N,1);
%         y_mtx = H_mtx * x_mtx;
%         y_r = reshape(y_mtx.',N*M,1);
% 
%         [est_bits,ite,x_est] = MPA_detector(N,M,M_mod,sigma_2(iesn0),data_grid,y_mtx,H_mtx,15);
% %         [est_bits,ite,x_est] = MPA_detector(N,M,M_mod,1e-4,data_grid,y_mtx,H_mtx,15);
%         errors = sum(xor(est_bits,trans_info_bit));
%         err_ber(iesn0) = errors/N_bits_perfram/N_fram + err_ber(iesn0);
        
        %% count NMSE
        NMSE_nume = 0;  
        NMSE_deno = 0;
        for kk = 0:(N-1)
            for ll = 1:M
                h_w = 0;
                h_w_hat = 0;
                for pp = 1:P
                    if (test==1)
                        h_w = h_w + Sampling_Function_v(N,kk+1,1,k_v_test(pp)) * Sampling_Function_t(M,ll,1,l_t_test(pp)) * h_test(pp) * exp(-2i*pi*k_v_test(pp)*l_t_test(pp)/N/M);
                    else
                        h_w = h_w + Sampling_Function_v(N,kk+1,1,k_v_init(pp)) * Sampling_Function_t(M,ll,1,l_t_init(pp)) * h_c_init(pp) * exp(-2i*pi*k_v_init(pp)*l_t_init(pp)/N/M);
                    end
                end
                for vv = 1:virtual_size
                    h_w_hat = h_w_hat + Sampling_Function_v(N,kk+1,1,k_v_hat(vv)) * Sampling_Function_t(M,ll,1,l_t_hat(vv)) * h_hat(vv) * exp(-2i*pi*k_v_hat(vv)*l_t_hat(vv)/N/M);
                end
                NMSE_nume = NMSE_nume + abs(h_w - h_w_hat).^2;
                NMSE_deno = NMSE_deno + abs(h_w)^2;
            end 
        end
        NMSE_count(iesn0) = NMSE_count(iesn0) + NMSE_nume / (NMSE_deno * N_fram);
        display(ifram,'ifram');
        display(iesn0, 'iesn0');
    end
end
NMSE_dB = 10 * log10(NMSE_count);

figure()

plot(SNR_dB,NMSE_dB,'-s','LineWidth',2,'MarkerSize',8)
grid on
xlabel('SNR(dB)')
ylabel('NMSE(dB)')

% figure;
% semilogy(SNR_dB, err_ber,'-*','LineWidth',2);
% grid on
% xlabel('SNR(dB)')
% ylabel('BER'); 