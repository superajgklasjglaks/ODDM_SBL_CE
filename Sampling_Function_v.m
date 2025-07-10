function w = Sampling_Function_v(N,k,k_in,k_v)
    pix = pi.*(k-k_in-k_v); 
    w = (sin(pix)./sin(pix./N).*exp(-1i*(N-1).*pix./N))/N;    % equation 12 & 13
    w(pix==0) = 1;
end

% function w = Sampling_Function_v(N,k,k_in,k_v)
%     w = 10*exp(1j*pi*k*sin(k_v))/N;
% end