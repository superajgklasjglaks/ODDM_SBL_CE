function A_beta0 = Gen_A_beta0(M_t,N_v,phi_trunc,alpha,beta0,mu_hbar,sigma_hbar,y_T)
    resid = y_T - phi_trunc * mu_hbar;
    res2 = norm(resid, 'fro')^2;
%     sigma_re = reshape(sigma_hbar_upgrate.',N_v*M_t,1);
%     bsum = 0;
%     for kk = 0:N_v - 1
%         for ll = 1:M_t
% %             bsum = bsum + 1 - sigma_hbar_upgrate(kk+1,ll)/alpha(kk * M_t + ll);
%             adder = sigma_hbar_upgrate(kk * M_t + ll,kk * M_t + ll)/(alpha_prev(kk * M_t + ll)+1e-16);
% %             if alpha_prev(kk * M_t + ll)>1e-5
%             if 1
%                 bsum = bsum + 1 - sigma_hbar_upgrate(kk * M_t + ll,kk * M_t + ll)/(alpha_prev(kk * M_t + ll)+1e-16);
%             else
%                 bsum = bsum + 1;
%             end
%         end
%     end
    flag = max(phi_trunc * mu_hbar);
    gamma1 = 1 - real(diag(sigma_hbar)) ./ (alpha + 1e-16);
    sumgam = sum(gamma1);
    A_beta0 = res2 + 1 / beta0 * sum(gamma1); % equation 51
end