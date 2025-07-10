function [phi_trunc,phi_t,phi_t_v,phi_t_t] = Gen_measurement_matrix(x_p,x_kp,x_lp,M,N,M_T,N_T,M_t,N_v,kv_bar,lt_bar,kappa,lorra)
%     phi_trunc = zeros(M_T*N_T,M_t*N_v);
    phi_t = zeros(M_T*N_T,M_t*N_v);
    phi_t_v = zeros(M_T*N_T,M_t*N_v);
    phi_t_t = zeros(M_T*N_T,M_t*N_v);
%     phi_t1 = zeros(M_T*N_T,M_t*N_v);
%     phi_t_t1 = zeros(M_T*N_T,M_t*N_v);
%     SF_v = zeros(M_T*N_T,M_t*N_v);
%     SF_t = zeros(M_T*N_T,M_t*N_v);
    for nn = 0:(N_v-1)
        for mm = 1:M_t
            for nnn = 0:(N_T-1)
                for mmm = 1:M_T
%                     aa = Sampling_Function(N,nnn,x_kp,kv_bar(nn * N_v + mm));
%                     bb = Sampling_Function(M,mmm,x_lp,lt_bar(nn * N_v + mm)); 
                    phi_t(nnn * M_T + mmm,nn * M_t + mm) = x_p * Sampling_Function_v(N,nnn,x_kp-1,kv_bar(nn * M_t + mm)) * Sampling_Function_t(M,mmm-1,x_lp-1,lt_bar(nn * M_t + mm)); 
                    phi_t_v(nnn * M_T + mmm,nn * M_t + mm) = x_p * wv_derivation(N,nnn,x_kp-1,kv_bar(nn * M_t + mm)) * Sampling_Function_t(M,mmm-1,x_lp-1,lt_bar(nn * M_t + mm));
                    phi_t_t(nnn * M_T + mmm,nn * M_t + mm) = x_p * Sampling_Function_v(N,nnn,x_kp-1,kv_bar(nn * M_t + mm)) * wt_derivation(M,mmm-1,x_lp-1,lt_bar(nn * M_t + mm));
%                     phi_t1(nnn * M_T + mmm,nn * M_t + mm) = x_p *  Sampling_Function_t(M,mmm-1,x_lp-1,lt_bar(nn * M_t + mm)); 
%                     phi_t_t1(nnn * M_T + mmm,nn * M_t + mm) = x_p * wt_derivation(M,mmm-1,x_lp-1,lt_bar(nn * M_t + mm)); 
%                     SF_v(nnn * M_T + mmm,nn * M_t + mm) = Sampling_Function_v(N,nnn,x_kp-1,kv_bar(nn * M_t + mm));
%                     SF_t(nnn * M_T + mmm,nn * M_t + mm) = Sampling_Function_t(M,mmm-1,x_lp-1,lt_bar(nn * M_t + mm));
%                     if(nnn * M_T + mmm==17)
%                         if(nn * M_t + mm==344)
%                             continue;
%                         end
%                     end
                end
            end
        end
    end
    phi_trunc = phi_t + phi_t_v * diag(kappa) + phi_t_t * diag(lorra);
%     lorra0 = lorra;
%     lorra0(50) = 0.25;
%     lorra0(51) = -0.25;
%     phi_test =  phi_t + phi_t_v * diag(kappa) + phi_t_t * diag(lorra0);
end

% function wvd = wv_derivation(N,k,k_p,k_v)
%     wvd = 0;
%     for n = 1:N-1
%        wvd = wvd + 2i*pi*(n-1)/N * exp(-2i*pi*(n-1)*(k-k_p-k_v)/N);   % equation 25
%     end
%     wvd = wvd/N;
% end
% 
% function wtd = wt_derivation(M,l,l_p,l_t)
%     wtd = 0;
%     for m = 1:M-1
%        wtd = wtd + 2i*pi*(m-1)/M * exp(2i*pi*(m-1)*(l-l_p-l_t)/M);   % equation 25
%     end
%     wtd = -wtd/M;
% end

%% x_p is x[x_kp,x_lp]
%{
kk = zeros(600,1);
kk(284) = 0;
ll = zeros(600,1);
ll(284) = 0.03;
phi_trunc = phi_t + phi_t_v * diag(kk) + phi_t_t * diag(ll);
%}