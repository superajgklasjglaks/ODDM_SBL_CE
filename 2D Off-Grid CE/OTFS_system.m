close all
clear all
rng(126)

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
r_v = 0.5;
r_t = 0.5;
P = 5;
% r_v = 1;
% r_t = 1;

%% parameters for test use
k_v_test = [1.06773328925145;0.0777528135573748;0.742234270223605;-0.137145209743403;-0.238441115602588];
l_t_test = [3.80345946429069;3.93830239631216;3.41769018357064;2.36436138031114;0.812785491218238];
h_test = [-0.228324281584917;-0.228178732316506;0.109874126055922;0.0373117918841507;-0.106173804509034];

% k_v_test = [1.06773328925145;0;0;0;0];
% l_t_test = [3.80345946429069;0;0;0;0];
% h_test = [-0.228324281584917;0;0;0;0];

% k_v_test = [0;0.0777528135573748;0;0;0];
% l_t_test = [0;3.93830239631216;0;0;0];
% h_test = [0;-0.228178732316506;0;0;0];

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
SNR_dB = 0:5:30;
% SNR_dB = 0:5:30;
SNR = 10.^(SNR_dB/10);
sigma_2 = 0.5 ./SNR;


%% Initializing simulation error count variables

N_fram = 100;
NMSE_count = zeros(length(SNR_dB),1);
% NMSE_dB = zeros(length(SNR_dB),1);
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
%     noise_re_init = mvnrnd(zeros(M*N,1),eye(M*N))';
%     noise_im_init = mvnrnd(zeros(M*N,1),eye(M*N))';
%         current_frame_number(iesn0)=ifram;
    for iesn0 = 1:length(SNR_dB)
        
        %% random input bits generation%%%%%
        trans_info_bit = randi([0,1],N_syms_perfram*M_bits,1);
        %%2D QAM symbols generation %%%%%%%%
        data=qammod(reshape(trans_info_bit,M_bits,N_syms_perfram), M_mod,'gray','InputType','bit');
        data = data./(eng_sqrt * 1);
        x = Generate_2D_data_grid_CE(N,M,data,data_grid);
        x(x_kp-floor(Kp/2):x_kp+floor(Kp/2),x_lp-floor(Lp/2):x_lp+floor(Lp/2))=sqrt(1000)/(Kp*Lp);
%         x(x_kp-floor(Kp/2):x_kp+floor(Kp/2),x_lp-floor(Lp/2):x_lp+floor(Lp/2))=1000/(Kp*Lp);
        
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
%         noise_gen_re = noise_re_init;
%         noise_gen_im = sigma_2(iesn0).*noise_im_init;
        noise_gen_re = mvnrnd(zeros(M*N,1),sigma_2(iesn0).*eye(M*N))';
        noise_gen_im = mvnrnd(zeros(M*N,1),sigma_2(iesn0).*eye(M*N))';
        noise_gen = noise_gen_re + 1i.* noise_gen_im;
%         noise_gen = noise_gen./10;
        h_test_hat = h_test .* exp(-2i * pi/M/N * (k_v_test.*l_t_test));
        h_c_init_hat = h_c_init .* exp(-2i * pi/M/N * (k_v_init.*l_t_init));
        if (test==1)
            r = phi_sys * h_test_hat;
%             r = phi_sys * h_test_hat + noise_gen;
        else
%             r = phi_sys * h_c_init_hat;
            r = phi_sys * h_c_init_hat + noise_gen;
        end
        y = reshape(r,M,N).';
        
        %% get truncation measument matrix%%%%%
%         y_trunc = zeros(2*k_max+Kp,l_max+Lp);
        y_trunc = y(floor(N/2)-floor(Kp/2)-k_max:floor(N/2)+floor(Kp/2)+k_max,floor(M/2)-floor(Lp/2):floor(M/2)+floor(Lp/2)+l_max);
        N_T = 2*k_max+Kp;
        M_T = l_max+Lp;
        y_T = reshape(y_trunc.',N_T*M_T,1);
%         [H_opt,k_v_opt,l_t_opt,virtual_size] = CE_2D_SBL(sqrt(sqrt(1000))/(Kp*Lp),k_max+Kp,sqrt(sqrt(1000))/(Kp*Lp),Lp,M,N,N_T,M_T,y_trunc,r_v,r_t,k_max,l_max);  % paper miss
        [H_opt,k_v_opt,l_t_opt,virtual_size] = CE_2D_SBL(sqrt(1000)/(2*Kp*Lp),k_max+Kp,2,Lp,M,N,N_T,M_T,y_trunc,r_v,r_t,k_max,l_max);  % paper miss
%         [H_opt,k_v_opt,l_t_opt,virtual_size] = CE_2D_SBL(sqrt(1000)/(Kp*Lp),k_max+Kp,sqrt(1000)/(Kp*Lp),Lp,M,N,N_T,M_T,y_trunc,r_v,r_t,k_max,l_max);  % paper miss
        
        %% count NMSE
        NMSE_nume = 0;  
        NMSE_deno = 0;
        for kk = 0:(N-1)
            for ll = 1:M
%                 h_w = 0;
%                 h_w_hat = 0;
%                 for pp = 1:P
%                     if (test==1)
%                         h_w = h_w + Sampling_Function_v(N,kk+1,1,k_v_test(pp)) * Sampling_Function_t(M,ll,1,l_t_test(pp)) * h_test(pp) * exp(-2i*pi*k_v_test(pp)*l_t_test(pp)/N/M);
%                     else
%                         h_w = h_w + Sampling_Function_v(N,kk+1,1,k_v_init(pp)) * Sampling_Function_t(M,ll,1,l_t_init(pp)) * h_c_init(pp) * exp(-2i*pi*k_v_init(pp)*l_t_init(pp)/N/M);
%                     end
%                 end
%                 for vv = 1:virtual_size
%                     h_w_hat = h_w_hat + Sampling_Function_v(N,kk+1,1,k_v_opt(vv)) * Sampling_Function_t(M,ll,1,l_t_opt(vv)) * H_opt(vv) * exp(-2i*pi*k_v_opt(vv)*l_t_opt(vv)/N/M);
%                 end
                if (test==1)
                    h_w = sum(Sampling_Function_v(N,kk+1,1,k_v_test) .* Sampling_Function_t(M,ll,1,l_t_test) .* h_test .* exp(-2i*pi.*k_v_test.*l_t_test/N/M));
                else
                    h_w = sum(Sampling_Function_v(N,kk+1,1,k_v_init) .* Sampling_Function_t(M,ll,1,l_t_init) .* h_c_init .* exp(-2i*pi.*k_v_init.*l_t_init/N/M));
                end
                h_w_hat = sum(Sampling_Function_v(N,kk+1,1,k_v_opt) .* Sampling_Function_t(M,ll,1,l_t_opt) .* H_opt .* exp(-2i*pi.*k_v_opt.*l_t_opt/N/M));
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

figure(1)

plot(SNR_dB,NMSE_dB,'-s','LineWidth',2,'MarkerSize',8)
grid on
xlabel('SNR(dB)')
ylabel('NMSE(dB)')
