function [H_opt, k_v_opt, l_t_opt] = hierarchical_SBL_refined_2D(x_v_p, x_kp, x_t_p, x_lp, M, N, N_T, M_T, Y_T, k_v_grid, l_t_grid, r_v, r_t, k_max, l_max, on_flag)
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