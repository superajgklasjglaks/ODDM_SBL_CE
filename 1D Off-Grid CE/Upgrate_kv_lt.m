function k_v_upgrate = Upgrate_kv_lt(k_v_prev,r_v)
%     label_low = find(k_v_prev<(-r_v)/2);
%     label_high = find(k_v_prev>(r_v)/2);
    k_v_upgrate = k_v_prev;
    k_v_upgrate(k_v_upgrate<(-r_v)/2) = (-r_v)/2;
    k_v_upgrate(k_v_upgrate>(r_v)/2) = (r_v)/2;
end
