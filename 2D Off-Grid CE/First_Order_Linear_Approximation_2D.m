function [kv_bar,lt_bar] = First_Order_Linear_Approximation_2D(N_v,M_t,v_max,r_v,r_t)
% N_v = 3;   % Doppler domain virtual sampling grid
% M_t = 3;   % delay domain virtual sampling grid
kv_bar = ones(N_v,1);
lt_bar = ones(M_t,1);
% k_v_init = [0:t_max/Nv:t_max];
% l_r_init = [-v_max:2*v_max/Mt:v_max];
for kk = 1:N_v
    kv_bar(kk) = (kk - 0.5) * r_v - v_max;
end
for ll = 1:M_t
    lt_bar(ll) = (ll - 0.5) * r_t;
end
% for kk = 0:N_v-1
%     for ll = 1:M_t
%         kv_bar(kk*M_t+ll) = (kk) * r_v - v_max;
%         lt_bar(kk*M_t+ll) = (ll - 1) * r_t;
%     end
% end
% 
% end

%%  First_Order_Linear_Approximation(5,3)