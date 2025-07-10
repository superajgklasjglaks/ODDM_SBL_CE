function [h_trd,k_v_trd,l_t_trd] = traditional_impulse(x_p,y_trunc,k_max,l_max,sigma_2_p)
    size_t = (2 * k_max + 1) * (l_max + 1);
    %H_bar = y_trunc*x_p / (x_p^2 + sigma_2_p*1e5);
    H_bar = y_trunc / x_p;
    H_bar(abs(H_bar)<(sigma_2_p * 3)) = 0;

    %% algorithm output
    h_trd = (reshape(H_bar.',1,size_t).');
    k_vv = repmat(-k_max:k_max,(l_max + 1),1);
    k_v_trd = (reshape(k_vv,1,size_t).');
    l_tt = repmat(0:l_max,(2 * k_max + 1),1);
    l_t_trd = (reshape(l_tt.',1,size_t).');
end