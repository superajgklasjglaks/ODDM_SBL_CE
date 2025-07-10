function [H_opt,k_v_opt,l_t_opt,virtual_size] = CE_2D_SBL(x_v_p,x_kp,x_t_p,x_lp,M,N,N_T,M_T,Y_T,r_v,r_t,k_max,l_max)
% kv_bar = ones(Nv * Mt,1);
% lt_bar = ones(Nv * Mt,1);
N_v = ceil(2 * k_max / r_v);   % Doppler domain virtual sampling grid
M_t = ceil(l_max / r_t);   % delay domain virtual sampling grid
virtual_size = N_v * M_t;
rho = 1e-2;
c = 1e-4;
d = 1e-4;
% zero_mu = zeros(virtual_size,1);
% alpha = gamrnd(1,row,virtual_size,1);   % equation 32
% delta_a = diag(alpha);
% h_bar = mvnrnd(zero_mu,delta_a);  % equation 31
% h_bar = h_bar';
% beta0 = gamrnd(1,d);   % equation 33
% z_T_bar = mvnrnd(zero_mu,1/beta0 * eye(virtual_size));
% z_T_bar = z_T_bar';
% k_v_prev = unifrnd(-r_v/2,r_v/2,virtual_size,1);
% l_t_prev = unifrnd(-r_t/2,r_t/2,virtual_size,1);

%% begin off-gird CE
% initialization
% y_T_normalized = y_T./x_p;
k_v = zeros(N_v,1);
l_t = zeros(M_t,1);
mu_dc = zeros(N_v,M_T);
A_beta0 = zeros(1,M_T);
% A_v = zeros(N_v,N_v);
% b_v = zeros(N_v,1);
H_bar = zeros(N_v,M_t);
l_t_bar = zeros(N_v,M_t);
H_opt = zeros(virtual_size,1);
k_v_opt = zeros(virtual_size,1);
l_t_opt = zeros(virtual_size,1);
[kv_bar,lt_bar] = First_Order_Linear_Approximation_2D(N_v,M_t,k_max,r_v,r_t);
[phi_trunc_L,phi_L,phi_L_v] = Gen_measurement_matrix_L(x_v_p,x_kp,N,N_T,N_v,kv_bar,k_v);
[~,phi_R,phi_R_t] = Gen_measurement_matrix_R(x_t_p,x_lp,M,M_T,M_t,lt_bar,l_t);

ksi = 1e-3;
% ksi = 1e-3;
% ksi = 1e-4;
Tmax = 2e2;
% test params
betas = zeros(Tmax,1);
toles = zeros(Tmax,1);
kv_test = zeros(Tmax,1);
% lt_test = zeros(Tmax,1);
mu_test = zeros(Tmax,1);
A_vs = zeros(Tmax,1);
b_vs = zeros(Tmax,1);
% b_v1s = zeros(Tmax,1);
% b_v2s = zeros(Tmax,1);
% A_T_ts = zeros(Tmax,1);
% b_T_ts = zeros(Tmax,1);

%% initialization
sigma2_bar = sum(sum(abs(Y_T).^2))/(100*N_T*M_T);
% P_hat = floor(M_T*N_T/log(virtual_size));
beta0 = 1/sigma2_bar;
phiLY = phi_trunc_L' * Y_T;
% alpha_v1 = zeros(N_v,1);
% for kpp = 1:N_v
%     for ll = 1:M_T
%         alpha_v1(kpp) = alpha_v1(kpp) + abs(phiLY(kpp,ll))/M_T;
%     end
% end
alpha_v = mean(abs(phiLY),2);
% da = alpha_v - alpha_v1;

%% begin algorithm 2
for t = 1:Tmax 
    phi_trunc_L = phi_L + phi_L_v * diag(k_v);
    alpha_v_prev = alpha_v;
    delta_av = diag(alpha_v);
    % algo 1 step 3
    % upgrate mu and Sigma
%     [mu_hbar,sigma_hbar] = Cal_upgrate_parameters(M_T,N_T,beta0,phi_trunc,y_T,delta_a);    % equation 40 & 41
    C = (1 / beta0) * eye(N_T) + phi_trunc_L * delta_av * phi_trunc_L';
    C_inv = inv(C);
    sigma_dc = delta_av - delta_av * phi_trunc_L' * C_inv * phi_trunc_L * delta_av;  % equation 41
    for lm = 1:M_T
        mu_dc(:,lm) = beta0 * sigma_dc * phi_trunc_L' * Y_T(:,lm);    % equation 40
    end
    
    % algo 1 step 4
    % count gamma1 for upgrating beta0
    gamma1 = 1 - real(diag(sigma_dc)) ./ (alpha_v + 1e-16);
    
    % upgrate alpha
    musq = sum(abs(mu_dc).^2, 2);
    alpha_temp = musq + M_T * diag(sigma_dc);
    alpha_v = -0.5 * M_T / rho + sqrt(0.25 * M_T^2 / rho^2 + alpha_temp /rho);  % equation 49
    
    % upgrate beta0
%     A_beta0 = Gen_A_beta0(M_t,N_v,phi_trunc,alpha,beta0,mu_hbar,sigma_hbar,y_T); % equation 51
    for lr = 1:M_T
        resid = Y_T(:,lr) - phi_trunc_L * mu_dc(:,lr);
        res2 = norm(resid)^2;
        A_beta0(lr) = res2 + 1 / beta0 * sum(gamma1); % equation 51
    end
%     ssum = sum(A_beta0);
    beta0 = (c - 1 + M_T * N_T)/(d + sum(A_beta0));
    
    % end condition
    tolerance = norm(alpha_v - alpha_v_prev)/norm(alpha_v_prev);
    if tolerance <= ksi
%         disp(t);
        break;
    end
    
    % upgrate matrix D

%     [~, idx] = sort(alpha_v, 'descend');
%     k_v_prev = k_v;
    
    if(tolerance<1000 * ksi)
%     
    %% upgrate k_v
    mmph = zeros(N_v,N_v);
    pLpv = phi_L_v' * phi_L_v;
    A_v = zeros(N_v,N_v);
    b_v = zeros(N_v,1);
    b_v = zeros(N_v,1);
    for l = 1:M_T
%         A_v = A_v + real(conj(pLpv) .* (mu_dc(:,l) * mu_dc(:,l)' + sigma_dc))/M_T; % equation 78
        A_v = A_v + real(pLpv .* (conj(mu_dc(:,l)) * mu_dc(:,l).' + sigma_dc.'))/M_T; % equation 78
%         b_v1 = diag(conj(mu_dc(:,l))) * phi_L_v' * (Y_T(:,l) - phi_L * mu_dc(:,l));
%         b_v = b_v + real(b_v1) / M_T;    % equation 55
        b_v1 = diag(mu_dc(:,l)) * phi_L_v.' * conj(Y_T(:,l)); 
        b_v2 = diag((mu_dc(:,l) * mu_dc(:,l)' + sigma_dc) * phi_L' * phi_L_v);
        b_v = b_v + real(b_v1 - b_v2) / M_T;    % equation 55
%         b_v3 = diag(mu_dc(:,l)) * phi_L_v.' * conj(Y_T(:,l)); 
%         b_v4 = diag((mu_dc(:,l) * mu_dc(:,l)' + sigma_dc) * phi_L' * phi_L_v);
%         b_vv = b_vv + real(b_v3 - b_v4) / M_T;    % equation 55
    end
%     b_v2 = diag(phi_L_v' * phi_L * sigma_dc);
%     b_v = b_v - real(b_v2);
%     bb = b_v - b_vv;
    
    temp = A_v \ b_v;
    k_v_upgrate = zeros(N_v,1);
    flag = ((max(svd(A_v)) / min(svd(A_v))) < 1e5) ;
    if (flag == 0) %|| any(abs(temp) > r_v/2)
        for pp = 1:N_v
            temp_k_v = k_v;
            temp_k_v(pp) = 0;
            if A_v(pp,pp) == 0
                k_v_upgrate(pp) = 0;
                continue;
            else
                k_v_upgrate(pp) = (b_v(pp) - A_v(pp,:) * temp_k_v)/A_v(pp,pp);   % equation 57
            end
        end
    else
        k_v_upgrate = temp;
    end
    
    k_v_upgrate(k_v_upgrate<(-r_v)/2) = (-r_v)/2;
    k_v_upgrate(k_v_upgrate>(r_v)/2) = (r_v)/2;
    k_v = k_v_upgrate;
    
    [~,pos_t] = max(alpha_v);
    
    if(t==100)
        zero = 0;
    end
    
    % monitor variables recording
%     lt_test(t) = l_t(pos_t);
    betas(t) = beta0;
    kv_test(t) = k_v(pos_t);
    mu_test(t) = mu_dc(pos_t);
    A_vs(t) = A_v(pos_t,pos_t);
    b_vs(t) = b_v(pos_t);
%     b_v1s(t) = b_v1(pos_t);
%     b_v2s(t) = b_v2(pos_t);
%     A_T_ts(t) = A_T_t(1,1);
%     b_T_ts(t) = b_T_t(1);
    
    end
    
    toles(t) = tolerance;
    % upgrate
%     mu_hbar = mu_hbar;
%     sigma_hbar = sigma_hbar; 
%     [phi_trunc,phi_t_v,phi_t_t] = Gen_measurement_matrix(x_p,x_kp,x_lp,M,N,M_T,N_T,M_t,N_v,kv_bar,lt_bar,k_v_prev,l_t_prev);
    
end

%% step 1 output
D_hat = mu_dc.';
k_v_hat = kv_bar + k_v;

%% step 2
for stp = 1:N_v
    [h_hat,l_t_hat] = CE_Algo1(M_T,D_hat(:,stp),r_t,M_t,phi_R,phi_R_t,lt_bar);
    H_bar(stp,:) = (h_hat.');
    l_t_bar(stp,:) = (l_t_hat.');
end

%% algorithm output
H_opt = (reshape(H_bar.',1,virtual_size).');
k_vv = repmat(k_v_hat,1,M_t);
k_v_opt = (reshape(k_vv.',1,virtual_size).');
l_t_opt = (reshape(l_t_bar.',1,virtual_size).');
end

%% CE_1D_SBL(3,3,ones(9,1),eye(9),eye(9),eye(9));