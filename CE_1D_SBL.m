function [h_hat,k_v_hat,l_t_hat,virtual_size, phi_trunc, delta_a] = CE_1D_SBL(x_p,x_kp,x_lp,M,N,N_T,M_T,y_T,r_v,r_t,k_max,l_max,on_flag)
N_v = ceil(2 * k_max / r_v);   % Doppler domain virtual sampling grid
M_t = ceil(l_max / r_t);   % delay domain virtual sampling grid
virtual_size = N_v * M_t;
rho = 1e-2;
c = 1e-4;
d = 1e-4;

%% begin off-gird CE
k_v = zeros(virtual_size,1);
l_t = zeros(virtual_size,1);
[kv_bar,lt_bar] = First_Order_Linear_Approximation(N_v,M_t,k_max,r_v,r_t);
[phi_trunc,phi_t,phi_t_v,phi_t_t] = Gen_measurement_matrix(x_p,x_kp,x_lp,M,N,M_T,N_T,M_t,N_v,kv_bar,lt_bar,k_v,l_t);

ksi = 1e-3;
Tmax = 5e2;

sigma2_bar = sum(abs(y_T).^2)/(100*N_T*M_T);
P_hat = floor(M_T*N_T/log(virtual_size));
alpha = abs(phi_trunc.' * y_T);
beta0 = 1/sigma2_bar;

% repeat
for t = 1:Tmax 
    phi_trunc = phi_t + phi_t_v * diag(k_v) + phi_t_t * diag(l_t);
    alpha_prev = alpha;
    delta_a = diag(alpha);
    % algo 1 step 3
    % upgrate mu and Sigma
    C = (1 / beta0) * eye(M_T*N_T) + phi_trunc * delta_a * phi_trunc';
    C_inv = inv(C);
    sigma_hbar = delta_a - delta_a * phi_trunc' * C_inv * phi_trunc * delta_a;  % equation 41
    mu_hbar = beta0 * sigma_hbar * phi_trunc' * y_T;    % equation 40    
    % algo 1 step 4
    % count gamma1 for upgrating beta0
    gamma1 = 1 - real(diag(sigma_hbar)) ./ (alpha + 1e-16);
    
    % upgrate alpha
    alpha_temp = abs(mu_hbar).^2 + real(diag(sigma_hbar));
    alpha = -0.5/rho + sqrt(0.25 / rho^2 + alpha_temp /rho);  % equation 49
    
    % upgrate beta0
    resid = y_T - phi_trunc * mu_hbar;
    res2 = norm(resid, 'fro')^2;
    A_beta0 = res2 + 1 / beta0 * sum(gamma1); % equation 51
    beta0 = (c - 1 + M_T * N_T)/(d + A_beta0);
    
    % end condition
    tolerance = norm(alpha - alpha_prev)/norm(alpha_prev);
    if tolerance <= ksi
        break;
    end


    [~, idx] = sort(alpha, 'descend');
    k_v_prev = k_v;
    l_t_prev = l_t;
    mmph = (mu_hbar * mu_hbar' + sigma_hbar);
    
    if(tolerance<1000 * ksi)
    %% upgrate k_v
    idx = idx(1:P_hat);
    
    pHpv = phi_t_v' * phi_t_v;
    tptlv = (phi_t + phi_t_t * diag(l_t_prev));
    A_v = real(pHpv .* (conj(mu_hbar) * mu_hbar.' + sigma_hbar.')); % equation 54
    b_v1 = diag(mu_hbar) * phi_t_v.' * conj(y_T);
    b_v2 = diag(mmph * tptlv' *  phi_t_v);
    b_v = real(b_v1 - b_v2);    % equation 55
    
    A_T_v = A_v(idx,idx);
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
    if(on_flag==1)
        k_v = k_v_upgrate;
    end
    
    %% upgrate l_t
    
    pHpt = phi_t_t' * phi_t_t;
    tptlt = (phi_t + phi_t_v * diag(k_v_prev));
    A_t = real(pHpt .* (conj(mu_hbar) * mu_hbar.' + sigma_hbar.')); % equation 54
    b_t1 = diag(mu_hbar) * phi_t_t.' * conj(y_T);
    b_t2 = diag((mmph * tptlt') *  phi_t_t);
    b_t = real(b_t1 - b_t2);    % equation 55
    
    A_T_t = A_t(idx,idx);
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
    if(on_flag==1)
        l_t = l_t_upgrate;
    end
    
    end
    
    toles(t) = tolerance;
end

%% output
h_hat = mu_hbar;
k_v_hat = kv_bar + k_v;
l_t_hat = lt_bar + l_t;
end