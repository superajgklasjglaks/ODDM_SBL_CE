function [phi_trunc_L,phi_L,phi_L_v] = Gen_measurement_matrix_L(x_v_p,x_kp,N,N_T,N_v,kv_bar,k_v)
    phi_L = zeros(N_T,N_v);
    phi_L_v = zeros(N_T,N_v);
    for nn = 1:N_v
        for nnn = 1:N_T
            phi_L(nnn,nn) = x_v_p * Sampling_Function_v(N,nnn,x_kp,kv_bar(nn)) ; 
            phi_L_v(nnn,nn) = x_v_p * wv_derivation(N,nnn,x_kp,kv_bar(nn));
        end
    end
    phi_trunc_L = phi_L + phi_L_v * diag(k_v);
end

%% x_p is x[x_kp,x_lp]
%{
kk = zeros(600,1);
kk(284) = 0;
ll = zeros(600,1);
ll(284) = 0.03;
phi_trunc = phi_t + phi_t_v * diag(kk) + phi_t_t * diag(ll);
%}