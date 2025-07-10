function [h_hat,l_t_hat] = CE_Algo1(M_T,y_T,r_t,M_t,phi_R,phi_R_t,lt_bar,on_flag)
rho = 1e-2;
c = 1e-4;
d = 1e-4;

%% begin off-gird CE
% initialization
l_t = zeros(M_t,1);

ksi = 1e-3;
% ksi = 1e-3;
% ksi = 1e-4;
Tmax = 2e2;
% test params
betas = zeros(Tmax,1);
toles = zeros(Tmax,1);
% kv_test = zeros(Tmax,1);
% lt_test = zeros(Tmax,1);
% mu_test = zeros(Tmax,1);
% A_T_vs = zeros(Tmax,1);
% b_T_vs = zeros(Tmax,1);
% b_v1s = zeros(Tmax,1);
% b_v2s = zeros(Tmax,1);
% A_T_ts = zeros(Tmax,1);
% b_T_ts = zeros(Tmax,1);

sigma2_bar = sum(abs(y_T).^2)/(100*M_T);
P_hat = M_t;
phi_trunc_R = phi_R + phi_R_t * diag(l_t);
alpha = abs(phi_trunc_R.' * y_T);
% delta_a = diag(alpha);
beta0 = 1/sigma2_bar;


% repeat
for t = 1:Tmax 
    phi_trunc_R = phi_R + phi_R_t * diag(l_t);
    alpha_prev = alpha;
    delta_a = diag(alpha);
    % algo 1 step 3
    % upgrate mu and Sigma
    C = (1 / beta0) * eye(M_T) + phi_trunc_R * delta_a * phi_trunc_R';
    C_inv = inv(C);
    sigma_hbar = delta_a - delta_a * phi_trunc_R' * C_inv * phi_trunc_R * delta_a;  % equation 41
    mu_hbar = beta0 * sigma_hbar * phi_trunc_R' * y_T;    % equation 40
    
    % algo 1 step 4
    % count gamma1 for upgrating beta0
    gamma1 = 1 - real(diag(sigma_hbar)) ./ (alpha + 1e-16);
    
    % upgrate alpha
%     alpha_upgrate = (sqrt(1 + 4 * row .* abs(mu_hbar).^2 + diag(sigma_hbar))-1)./(2*row);  % equation 49
    alpha_temp = abs(mu_hbar).^2 + diag(sigma_hbar);
    alpha = -0.5/rho + sqrt(0.25 / rho^2 + alpha_temp /rho);  % equation 49
    
    % upgrate beta0
%     A_beta0 = Gen_A_beta0(M_t,N_v,phi_trunc,alpha,beta0,mu_hbar,sigma_hbar,y_T); % equation 51
    resid = y_T - phi_trunc_R * mu_hbar;
    res2 = norm(resid, 'fro')^2;
    A_beta0 = res2 + 1 / beta0 * sum(gamma1); % equation 51
    beta0 = (c - 1 + M_T)/(d + A_beta0);
    
    % end condition
    tolerance = norm(alpha - alpha_prev)/norm(alpha_prev);
    if tolerance <= ksi
%         disp(t);
        break;
    end
    
    % upgrate k_v and l_t
%     [A_v,b_v] = Gen_Av_bv(phi_trunc,phi_t_v,phi_t_t,l_t,mu_hbar,sigma_hbar,y_T); % equation 54 & 55
%     [A_t,b_t] = Gen_Av_bv(phi_trunc,phi_t_t,phi_t_v,k_v,mu_hbar,sigma_hbar,y_T); % equation 60 & 61
%     k_v = Upgrate_k_T_v(P_hat,A_v,b_v,alpha,k_v,r_v);   % equation 58
%     l_t = Upgrate_k_T_v(P_hat,A_t,b_t,alpha,l_t,r_t);   % equation 64


%     [~, idx] = sort(alpha, 'descend');
    l_t_prev = l_t;
    mmph = (mu_hbar * mu_hbar' + sigma_hbar);
    
    if(tolerance<1000 * ksi)    
    %% upgrate l_t
    pHpt = phi_R_t' * phi_R_t;
%     A_t = real(pHpt .* (conj(mu_hbar) * mu_hbar.' + conj(sigma_hbar))); % equation 54
    A_t = real(pHpt .* (conj(mu_hbar) * mu_hbar.' + sigma_hbar.')); % equation 54
    b_t1 = diag(mu_hbar) * phi_R_t.' * conj(y_T);
    b_t2 = diag((mmph * phi_R') *  phi_R_t);
%     A_t = real(conj(pHpt) .* (mu_hbar * mu_hbar' + sigma_hbar));
%     b_t1 = diag(conj(mu_hbar)) * phi_R_t' * (y_T - phi_R * mu_hbar);
%     b_t2 = diag(phi_R_t' * phi_R * sigma_hbar);
    b_t = real(b_t1 - b_t2);    % equation 55
    
    temp = A_t \ b_t;
    l_t_upgrate = zeros(M_t,1);
%     l_t_upgrate = l_t_prev;
    flag = ((max(svd(A_t)) / min(svd(A_t))) < 1e6);
    if (flag == 0) %|| any(abs(l_T_t) > r_t/2)
        for pp = 1:P_hat
            temp_l_t = l_t;
            temp_l_t(pp) = 0;
            if A_t(pp,pp) == 0
                l_t_upgrate(pp) = 0;
                continue;
            else
                l_t_upgrate(pp) = (b_t(pp) - A_t(pp,:) * temp_l_t)/A_t(pp,pp);   % equation 57
            end
        end
    else
        l_t_upgrate = temp;
    end
    l_t_upgrate(l_t_upgrate<(-r_t)/2) = (-r_t)/2;
    l_t_upgrate(l_t_upgrate>(r_t)/2) = (r_t)/2;
    if(on_flag==1)
        l_t = l_t_upgrate;
    end
    
%     [~,pos_t] = max(alpha);
    
    % monitor variables recording
%     betas(t) = A_beta0;
%     lt_test(t) = l_t(pos_t);
%     mu_test(t) = mu_hbar(pos_t);
%     A_T_ts(t) = A_T_t(1,1);
%     b_T_ts(t) = b_T_t(1);
    
    end
    
    toles(t) = tolerance;
    
end

%% output
h_hat = mu_hbar;
l_t_hat = lt_bar + l_t;

end

%% CE_1D_SBL(3,3,ones(9,1),eye(9),eye(9),eye(9));