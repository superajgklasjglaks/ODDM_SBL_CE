function wvd = wv_derivation(N,k,k_p,k_v)
    wvd = 0;
    for n = 1:N
       wvd = wvd + 2i*pi*n/N * exp(-2i*pi*n.*(k-k_p-k_v)/N);   % equation 25
    end
    wvd = wvd/N;
end