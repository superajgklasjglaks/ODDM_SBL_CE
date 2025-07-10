function [A_v,b_v] = Gen_Av_bv(phi_trunc,phi_t_v,phi_t_t,lorra,mu_hbar_upgrate,sigma_hbar_upgrate,y_T)
%     A_v = real(conj(phi_t_v' * phi_t_v) .* (mu_hbar_upgrate * mu_hbar_upgrate' + sigma_hbar_upgrate)); % equation 54
    A_v = real((phi_t_v' * phi_t_v) .* (conj(mu_hbar_upgrate) * mu_hbar_upgrate.' + sigma_hbar_upgrate.')); % equation 54
    b_v = diag(mu_hbar_upgrate) * phi_t_v.' * conj(y_T);
    b_v = b_v - diag((mu_hbar_upgrate * mu_hbar_upgrate' + sigma_hbar_upgrate)*(phi_trunc + phi_t_t * diag(lorra))' *  phi_t_v);
    b_v = real(b_v);    % equation 55
end

%%
%{
max(max(abs(A_v)))
max(max(abs(phi_t_v' * phi_t_v)))
max(max(abs((conj(mu_hbar_upgrate) * mu_hbar_upgrate.' + sigma_hbar_upgrate.'))))
%}