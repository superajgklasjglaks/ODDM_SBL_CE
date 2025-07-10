function k_v_upgrate = Upgrate_k_v_o(virtual_size,P_hat,idx,phi_t,phi_t_v,phi_t_t,lorra,mu_hbar,sigma_hbar,y_T,k_v,r_v)
    idx = idx(1:P_hat);
    
    pHp = phi_t_v' * phi_t_v;
    A_v = real(pHp .* (conj(mu_hbar) * mu_hbar.' + sigma_hbar.')); % equation 54
    b_v = diag(mu_hbar) * phi_t_v.' * conj(y_T);
    mmph = (mu_hbar * mu_hbar' + sigma_hbar);
    mmph_m = max(max(mmph));
    tptl = (phi_t + phi_t_t * diag(lorra));
    tptl_m = max(max(tptl));
    b_v = b_v - diag((mmph * tptl') *  phi_t_v);
    b_v = real(b_v);    % equation 55
    
    A_T_v = A_v(idx,idx);
    A_T_v = diag(diag(A_T_v));
    b_T_v = b_v(idx);
    
    temp = k_v;
    k_v_upgrate = zeros(virtual_size,1);
    k_v_upgrate(idx) = temp(idx);
    flag = ((max(svd(A_T_v)) / min(svd(A_T_v))) < 1e6);
    if (flag == 0)
        for n = 1:length(k_v)
            if A_T_v(n,n) == 0
                k_v_upgrate(idx(n)) = 0;
                continue;
            else
                for pp = 1:P_hat
                     k_v(pp) = (b_T_v(pp) - A_T_v(:,pp)' * k_v +  A_T_v(pp,pp) * k_v(pp))/A_T_v(pp,pp);   % equation 57
                end
            end
        end
    else
        k_v_T = A_T_v \ b_T_v;
        k_v_upgrate(idx) = A_T_v \ b_T_v;
    end
    k_v_upgrate(k_v_upgrate<(-r_v)/2) = (-r_v)/2;
    k_v_upgrate(k_v_upgrate>(r_v)/2) = (r_v)/2;
end
