function [h_hat, k_v_hat, l_t_hat] = hierarchical_SBL_refined_1D(x_p, x_kp, x_lp, M, N, N_T, M_T, y_T, k_v_grid, l_t_grid, on_flag)
    % Custom SBL implementation for predefined grid points
    virtual_size = length(k_v_grid);
    rho = 1e-2;
    c = 1e-4;
    d = 1e-4;
    
    % Create measurement matrix for custom grid
    phi_trunc = zeros(M_T*N_T, virtual_size);
    for i = 1:virtual_size
        for nnn = 0:(N_T-1)
            for mmm = 1:M_T
                phi_trunc(nnn * M_T + mmm, i) = x_p * Sampling_Function_v(N,nnn,x_kp-1,k_v_grid(i)) * Sampling_Function_t(M,mmm-1,x_lp-1,l_t_grid(i));
            end
        end
    end
    
    ksi = 1e-3;
    Tmax = 5e2;
    
    sigma2_bar = sum(abs(y_T).^2)/(100*N_T*M_T);
    P_hat = floor(M_T*N_T/log(virtual_size));
    alpha = abs(phi_trunc.' * y_T);
    beta0 = 1/sigma2_bar;
    
    % SBL iteration
    for t = 1:Tmax
        alpha_prev = alpha;
        delta_a = diag(alpha);
        
        % Update mu and Sigma
        C = (1 / beta0) * eye(M_T*N_T) + phi_trunc * delta_a * phi_trunc';
        C_inv = inv(C);
        sigma_hbar = delta_a - delta_a * phi_trunc' * C_inv * phi_trunc * delta_a;
        mu_hbar = beta0 * sigma_hbar * phi_trunc' * y_T;
        
        % Update gamma1
        gamma1 = 1 - real(diag(sigma_hbar)) ./ (alpha + 1e-16);
        
        % Update alpha
        alpha_temp = abs(mu_hbar).^2 + real(diag(sigma_hbar));
        alpha = -0.5/rho + sqrt(0.25 / rho^2 + alpha_temp /rho);
        
        % Update beta0
        resid = y_T - phi_trunc * mu_hbar;
        res2 = norm(resid, 'fro')^2;
        A_beta0 = res2 + 1 / beta0 * sum(gamma1);
        beta0 = (c - 1 + M_T * N_T)/(d + A_beta0);
        
        % Check convergence
        tolerance = norm(alpha - alpha_prev)/norm(alpha_prev);
        if tolerance <= ksi
            break;
        end
    end
    
    % Output
    h_hat = mu_hbar;
    k_v_hat = k_v_grid;
    l_t_hat = l_t_grid;
end