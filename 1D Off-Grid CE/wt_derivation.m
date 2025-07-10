function wtd = wt_derivation(M,l,l_p,l_t)
    adder = 0;
    for m = 1:M
       adder = adder + 2i*pi*m/M * exp(2i*pi*m.*(l-l_p-l_t)/M);   % equation 25
    end
    wtd = -adder/M;
end
