function [phi_trunc,phi_t,phi_t_v,phi_t_t] = Gen_measurement_matrix(x_p,x_kp,x_lp,M,N,M_T,N_T,M_t,N_v,kv_bar,lt_bar,kappa,lorra)
    phi_t = zeros(M_T*N_T,M_t*N_v);
    phi_t_v = zeros(M_T*N_T,M_t*N_v);
    phi_t_t = zeros(M_T*N_T,M_t*N_v);
    for nn = 0:(N_v-1)
        for mm = 1:M_t
            for nnn = 0:(N_T-1)
                for mmm = 1:M_T
                    phi_t(nnn * M_T + mmm,nn * M_t + mm) = x_p * Sampling_Function_v(N,nnn,x_kp-1,kv_bar(nn * M_t + mm)) * Sampling_Function_t(M,mmm-1,x_lp-1,lt_bar(nn * M_t + mm)); 
                    phi_t_v(nnn * M_T + mmm,nn * M_t + mm) = x_p * wv_derivation(N,nnn,x_kp-1,kv_bar(nn * M_t + mm)) * Sampling_Function_t(M,mmm-1,x_lp-1,lt_bar(nn * M_t + mm));
                    phi_t_t(nnn * M_T + mmm,nn * M_t + mm) = x_p * Sampling_Function_v(N,nnn,x_kp-1,kv_bar(nn * M_t + mm)) * wt_derivation(M,mmm-1,x_lp-1,lt_bar(nn * M_t + mm));
                end
            end
        end
    end
    phi_trunc = phi_t + phi_t_v * diag(kappa) + phi_t_t * diag(lorra);
end
