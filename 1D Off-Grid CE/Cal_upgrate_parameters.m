function [mu_hbar_upgrate,sigma_hbar_upgrate] = Cal_upgrate_parameters(M_T,N_T,beta0,phi_trunc,y_T,delta_a)
    C = 1 / beta0 * eye(M_T*N_T) + phi_trunc * delta_a * phi_trunc';
    C_inv = inv(C);
    sigma_hbar_upgrate = delta_a - delta_a * phi_trunc' * C_inv * phi_trunc * delta_a;  % equation 41
    mu_hbar_upgrate = beta0 * sigma_hbar_upgrate * phi_trunc' * y_T;    % equation 40
    delta_upg = delta_a * phi_trunc' * inv((1/beta0).* eye(M_T*N_T) + phi_trunc * delta_a * phi_trunc') * phi_trunc * delta_a;
%     [sigma_max1,pos1] = max(delta_upg);
%     [sigma_max2,pos2] = max(sigma_max1);
    sigma_hbar_upgrate_1 = inv(beta0.* phi_trunc' * phi_trunc + inv(delta_a));
    diff = sigma_hbar_upgrate - sigma_hbar_upgrate_1;
    abs_diff = max(max(abs(diff)));
%     mu_hbar_upgrate(mu_hbar_upgrate<1e-5) = 0;
%     sigma_hbar_upgrate(abs(sigma_hbar_upgrate).^2<1e-5) = 0;
end
