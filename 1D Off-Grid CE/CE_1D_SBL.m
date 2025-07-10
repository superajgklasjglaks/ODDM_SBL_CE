function [h_hat,k_v_hat,l_t_hat,virtual_size] = CE_1D_SBL(x_p,x_kp,x_lp,M,N,N_T,M_T,y_T,r_v,r_t,k_max,l_max)
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
k_v = zeros(virtual_size,1);
l_t = zeros(virtual_size,1);
% idx = [];
[kv_bar,lt_bar] = First_Order_Linear_Approximation(N_v,M_t,k_max,r_v,r_t);
[phi_trunc,phi_t,phi_t_v,phi_t_t] = Gen_measurement_matrix(x_p,x_kp,x_lp,M,N,M_T,N_T,M_t,N_v,kv_bar,lt_bar,k_v,l_t);

ksi = 1e-3;
% ksi = 1e-3;
% ksi = 1e-4;
Tmax = 5e2;
% test params
% betas = zeros(Tmax,1);
% toles = zeros(Tmax,1);
% kv_test = zeros(Tmax,1);
% lt_test = zeros(Tmax,1);
% mu_test = zeros(Tmax,1);
% A_T_vs = zeros(Tmax,1);
% b_T_vs = zeros(Tmax,1);
% b_v1s = zeros(Tmax,1);
% b_v2s = zeros(Tmax,1);
% A_T_ts = zeros(Tmax,1);
% b_T_ts = zeros(Tmax,1);

sigma2_bar = sum(abs(y_T).^2)/(100*N_T*M_T);
P_hat = floor(M_T*N_T/log(virtual_size));
% P_hat = 1;
alpha = abs(phi_trunc.' * y_T);
% delta_a = diag(alpha);
beta0 = 1/sigma2_bar;

% sigma_hbar = delta_a;
% sigma_hbar = sigma_hbar - delta_a * phi_trunc' * (((1/beta0).* eye(M_T*N_T) + phi_trunc * delta_a * phi_trunc') \ phi_trunc) * delta_a;  % equation 41;
% delta_upg = delta_a * phi_trunc' * (((1/beta0).* eye(M_T*N_T) + phi_trunc * delta_a * phi_trunc') \ phi_trunc) * delta_a;
% mu_hbar = beta0 .* sigma_hbar * phi_trunc' * y_T;    % equation 40

% repeat
for t = 1:Tmax 
%     phi_trunc(:,idx) = phi_t(:,idx) + phi_t_v(:,idx) * diag(k_v(idx)) + phi_t_t(:,idx) * diag(l_t(idx));
    phi_trunc = phi_t + phi_t_v * diag(k_v) + phi_t_t * diag(l_t);
    alpha_prev = alpha;
    delta_a = diag(alpha);
    % algo 1 step 3
    % upgrate mu and Sigma
%     [mu_hbar,sigma_hbar] = Cal_upgrate_parameters(M_T,N_T,beta0,phi_trunc,y_T,delta_a);    % equation 40 & 41
    C = (1 / beta0) * eye(M_T*N_T) + phi_trunc * delta_a * phi_trunc';
    C_inv = inv(C);
    sigma_hbar = delta_a - delta_a * phi_trunc' * C_inv * phi_trunc * delta_a;  % equation 41
    mu_hbar = beta0 * sigma_hbar * phi_trunc' * y_T;    % equation 40
%     flag_mu_sigma = zeros(virtual_size,1);
    
    % algo 1 step 4
    % count gamma1 for upgrating beta0
    gamma1 = 1 - real(diag(sigma_hbar)) ./ (alpha + 1e-16);
    
    % upgrate alpha
%     alpha_upgrate = (sqrt(1 + 4 * row .* abs(mu_hbar).^2 + diag(sigma_hbar))-1)./(2*row);  % equation 49
    alpha_temp = abs(mu_hbar).^2 + real(diag(sigma_hbar));
    alpha = -0.5/rho + sqrt(0.25 / rho^2 + alpha_temp /rho);  % equation 49
    
    % upgrate beta0
%     A_beta0 = Gen_A_beta0(M_t,N_v,phi_trunc,alpha,beta0,mu_hbar,sigma_hbar,y_T); % equation 51
    resid = y_T - phi_trunc * mu_hbar;
    res2 = norm(resid, 'fro')^2;
    A_beta0 = res2 + 1 / beta0 * sum(gamma1); % equation 51
    beta0 = (c - 1 + M_T * N_T)/(d + A_beta0);
    
    % end condition
    tolerance = norm(alpha - alpha_prev)/norm(alpha_prev);
    if tolerance <= ksi
        disp(t);
        break;
    end
    
    % upgrate k_v and l_t
%     [A_v,b_v] = Gen_Av_bv(phi_trunc,phi_t_v,phi_t_t,l_t,mu_hbar,sigma_hbar,y_T); % equation 54 & 55
%     [A_t,b_t] = Gen_Av_bv(phi_trunc,phi_t_t,phi_t_v,k_v,mu_hbar,sigma_hbar,y_T); % equation 60 & 61
%     k_v = Upgrate_k_T_v(P_hat,A_v,b_v,alpha,k_v,r_v);   % equation 58
%     l_t = Upgrate_k_T_v(P_hat,A_t,b_t,alpha,l_t,r_t);   % equation 64


    [~, idx] = sort(alpha, 'descend');
    k_v_prev = k_v;
    l_t_prev = l_t;
    mmph = (mu_hbar * mu_hbar' + sigma_hbar);
    
    if(tolerance<1000 * ksi)
%     
    %% upgrate k_v
    idx = idx(1:P_hat);
    
    pHpv = phi_t_v' * phi_t_v;
    tptlv = (phi_t + phi_t_t * diag(l_t_prev));
%     A_v = real(pHpv .* (conj(mu_hbar) * mu_hbar.' + conj(sigma_hbar))); % equation 54
    A_v = real(pHpv .* (conj(mu_hbar) * mu_hbar.' + sigma_hbar.')); % equation 54
    b_v1 = diag(mu_hbar) * phi_t_v.' * conj(y_T);
    b_v2 = diag(mmph * tptlv' *  phi_t_v);
%     A_v = real(conj(pHpv) .* (mu_hbar * mu_hbar' + sigma_hbar));
%     b_v1 = real(diag(conj(mu_hbar)) * phi_t_v' * (y_T - tptlv * mu_hbar));
% %     deltay = y_T - tptlv * mu_hbar; 
%     b_v2 = real(diag(phi_t_v' * tptlv * sigma_hbar));
    b_v = real(b_v1 - b_v2);    % equation 55
    
    A_T_v = A_v(idx,idx);
%     A_T_v = diag(diag(A_T_v));
    b_T_v = b_v(idx);
    
    temp = k_v;
    k_T_v = A_T_v \ b_T_v;
    k_v_upgrate = zeros(virtual_size,1);
    k_v_upgrate(idx) = temp(idx);
    flag = ((max(svd(A_T_v)) / min(svd(A_T_v))) < 1e6) ;
    if (flag == 0) %|| any(abs(k_T_v) > r_v/2)
        for pp = 1:P_hat
            temp_k_v = k_v(idx);
            temp_k_v(pp) = 0;
            if A_T_v(pp,pp) == 0
                k_v_upgrate(idx(pp)) = 0;
                continue;
            else
                k_v_upgrate(idx(pp)) = (b_T_v(pp) - A_T_v(pp,:) * temp_k_v)/A_T_v(pp,pp);   % equation 57
            end
        end
    else
        k_v_upgrate = zeros(virtual_size,1);
        k_v_upgrate(idx) = k_T_v;
    end
    
    k_v_upgrate(k_v_upgrate<(-r_v)/2) = (-r_v)/2;
    k_v_upgrate(k_v_upgrate>(r_v)/2) = (r_v)/2;
    k_v = k_v_upgrate;
    
    %% upgrate l_t
    
    pHpt = phi_t_t' * phi_t_t;
    tptlt = (phi_t + phi_t_v * diag(k_v_prev));
%     A_t = real(pHpt .* (conj(mu_hbar) * mu_hbar.' + conj(sigma_hbar))); % equation 54
    A_t = real(pHpt .* (conj(mu_hbar) * mu_hbar.' + sigma_hbar.')); % equation 54
    b_t1 = diag(mu_hbar) * phi_t_t.' * conj(y_T);
    b_t2 = diag((mmph * tptlt') *  phi_t_t);
%     A_t = real(conj(pHpt) .* (mu_hbar * mu_hbar' + sigma_hbar));
%     b_t1 = real(diag(conj(mu_hbar)) * phi_t_t' * (y_T - tptlt * mu_hbar));
%     b_t2 = real(diag(phi_t_t' * tptlt * sigma_hbar));
    b_t = real(b_t1 - b_t2);    % equation 55
    
    A_T_t = A_t(idx,idx);
%     A_T_t = diag(diag(A_T_t));
    b_T_t = b_t(idx);
    
    temp = l_t;
    l_T_t = A_T_t \ b_T_t;
    l_t_upgrate = zeros(virtual_size,1);
    l_t_upgrate(idx) = temp(idx);
    flag = ((max(svd(A_T_t)) / min(svd(A_T_t))) < 1e6);
    if (flag == 0) %|| any(abs(l_T_t) > r_t/2)
        for pp = 1:P_hat
            temp_l_t = l_t(idx);
            temp_l_t(pp) = 0;
            if A_T_t(pp,pp) == 0
                l_t_upgrate(idx(pp)) = 0;
                continue;
            else
                l_t_upgrate(idx(pp)) = (b_T_t(pp) - A_T_t(pp,:) * temp_l_t)/A_T_t(pp,pp);   % equation 57
            end
        end
    else
        l_t_upgrate = zeros(virtual_size,1);
        l_t_upgrate(idx) = l_T_t;
    end
    l_t_upgrate(l_t_upgrate<(-r_t)/2) = (-r_t)/2;
    l_t_upgrate(l_t_upgrate>(r_t)/2) = (r_t)/2;
    
    l_t = l_t_upgrate;
    
%     [~,pos_t] = max(alpha);
%     
%     if(t==100)
%         zero = 0;
%     end
    
    % monitor variables recording
%     betas(t) = A_beta0;
%     lt_test(t) = l_t(pos_t);
%     kv_test(t) = k_v(pos_t);
%     mu_test(t) = mu_hbar(pos_t);
%     A_T_vs(t) = A_T_v(1,1);
%     b_T_vs(t) = b_T_v(1);
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

%% output
% figure;
% betas(1:10)=0;
% plot(1:Tmax,betas);
% plot(1:Tmax,kv_test);
% figure;
% plot(1:Tmax,lt_test);
% figure;
% plot(1:Tmax,mu_test);
h_hat = mu_hbar;
k_v_hat = kv_bar + k_v;
l_t_hat = lt_bar + l_t;
% h_avg = sum(abs(h_hat).^2);

% disp('CE_1D_SBL Finished');

end

%% CE_1D_SBL(3,3,ones(9,1),eye(9),eye(9),eye(9));