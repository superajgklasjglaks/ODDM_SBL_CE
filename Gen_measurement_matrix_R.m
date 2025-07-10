function [phi_trunc_R,phi_R,phi_R_t] = Gen_measurement_matrix_R(x_t_p,x_lp,M,M_T,M_t,lt_bar,l_t)
    phi_R = zeros(M_T,M_t);
    phi_R_t = zeros(M_T,M_t);
    for mm = 1:M_t
        for mmm = 1:M_T
            phi_R(mmm,mm) = x_t_p * Sampling_Function_t(M,mmm,x_lp,lt_bar(mm)) ; 
            phi_R_t(mmm,mm) = x_t_p * wt_derivation(M,mmm,x_lp,lt_bar(mm));
        end
    end
    phi_trunc_R = phi_R + phi_R_t * diag(l_t);
end