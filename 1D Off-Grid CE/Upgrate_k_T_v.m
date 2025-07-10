function k_v_upgrate = Upgrate_k_T_v(P_hat,A_v,b_v,alpha,k_v,r_v)
    A_T_v = zeros(P_hat);
    b_T_v = zeros(P_hat,1);
    k_T_v = zeros(P_hat,1);
    S_v = zeros(P_hat,1);
    k_v_upgrate = zeros(length(k_v),1);
    alpha_del = alpha;
    for pp = 1:P_hat
        [~,pos] = max(alpha_del);
        alpha_del(pos) = 0;
        A_T_v(pp,pp) = A_v(pos,pos);
        b_T_v(pp) = b_v(pos);
        k_T_v(pp) = k_v(pos);
        S_v(pp) = pos;
    end
    flag = ((max(svd(A_T_v)) / min(svd(A_T_v))) < 1e4);
    if(flag~=0)
        k_T_v = A_T_v\b_T_v;  % equation 56
    else
        for pp = 1:P_hat
            k_T_v(pp) = (b_T_v(pp) - A_T_v(:,pp)' * k_T_v +  A_T_v(pp,pp) * k_T_v(pp))/A_T_v(pp,pp);   % equation 57
        end
    end
    
    k_T_v_upgrate = Upgrate_kv_lt(k_T_v,r_v);  % equation 58
    
    for pp = 1:P_hat
        k_v_upgrate(S_v(pp)) = k_T_v_upgrate(pp);
    end
end
