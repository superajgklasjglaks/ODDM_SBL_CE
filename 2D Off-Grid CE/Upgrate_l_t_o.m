function l_t_upgrate = Upgrate_l_t_o(virtual_size,P_hat,idx,phi_t,phi_t_t,phi_t_v,k_v,mu_hbar,sigma_hbar,y_T,l_t,r_t)
    idx = idx(1:P_hat);
    
    pHpt = phi_t_t' * phi_t_t;
    A_t = real(pHpt .* (conj(mu_hbar) * mu_hbar.' + sigma_hbar.')); % equation 54
    b_t = diag(mu_hbar) * phi_t_t.' * conj(y_T);
    mmpht = (mu_hbar * mu_hbar' + sigma_hbar);
%     mu_target = mu_hbar_upgrate(304);
%     sigma_target = sigma_hbar_upgrate(304,304);
%     mmph_m = max(max(mmpht));
    tptlt = (phi_t + phi_t_v * diag(k_v));
%     tptl_m = max(max(tptl));
    b_t = b_t - diag((mmpht * tptlt') *  phi_t_t);
    b_t = real(b_t);    % equation 55
    
    A_T_t = A_t(idx,idx);
    A_T_t = diag(diag(A_T_t));
    b_T_t = b_t(idx);
    
    temp = l_t;
    l_t_upgrate = zeros(virtual_size,1);
    l_t_upgrate(idx) = temp(idx);
    flag = ((max(svd(A_T_t)) / min(svd(A_T_t))) < 1e6);
    if (flag == 0)
        for n = 1:length(l_t)
            if A_T_t(n,n) == 0
                l_t_upgrate(idx(n)) = 0;
                continue;
            else
                for pp = 1:P_hat
                     l_t(pp) = (b_T_t(pp) - A_T_t(:,pp)' * l_t +  A_T_t(pp,pp) * l_t(pp))/A_T_t(pp,pp);   % equation 57
                end
            end
        end
    else
        k_v_T = A_T_t \ b_T_t;
        l_t_upgrate(idx) = A_T_t \ b_T_t;
    end
    l_t_upgrate(l_t_upgrate<(-r_t)/2) = (-r_t)/2;
    l_t_upgrate(l_t_upgrate>(r_t)/2) = (r_t)/2;
end
