function k_v_upgrate = Upgrate_k_v(virtual_size,P_hat,idx,phi_trunc,phi_t_v,phi_t_t,lorra,mu_hbar_upgrate,sigma_hbar_upgrate,y_T,alpha,k_v,r_v)
    idx = idx(1:P_hat);
    
    pHp = phi_t_v' * phi_t_v;
    A_v = real(pHp(idx,idx) .* (conj(mu_hbar_upgrate(idx)) * mu_hbar_upgrate(idx).' + sigma_hbar_upgrate(idx,idx).')); % equation 54
    b_v = diag(mu_hbar_upgrate(idx)) * phi_t_v(:,idx).' * conj(y_T);
    b_v = b_v - diag((mu_hbar_upgrate(idx) * mu_hbar_upgrate(idx)' + sigma_hbar_upgrate(idx,idx))*(phi_trunc(:,idx) + phi_t_t(:,idx) * diag(lorra(idx)))' *  phi_t_v(:,idx));
    b_v = real(b_v);    % equation 55
    
    temp = k_v;
    k_v_upgrate = zeros(virtual_size,1);
    k_v_upgrate(idx) = temp(idx);
    flag = ((max(svd(A_v)) / min(svd(A_v))) < 1e6);
    if (flag == 0)
        for n = 1:length(k_v)
            if A_v(n,n) == 0
                k_v_upgrate(idx(n)) = 0;
                continue;
            else
                for pp = 1:P_hat
                     k_v(pp) = (b_v(pp) - A_v(:,pp)' * k_v +  A_v(pp,pp) * k_v(pp))/A_v(pp,pp);   % equation 57
                end
            end
        end
    else
        k_v_upgrate(idx) = A_v \ b_v;
    end
    k_v_upgrate(k_v_upgrate<(-r_v)/2) = (-r_v)/2;
    k_v_upgrate(k_v_upgrate>(r_v)/2) = (r_v)/2;
end
