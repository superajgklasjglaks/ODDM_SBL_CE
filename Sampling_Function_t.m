% function w = Sampling_Function_t(M,l,l_in,l_t)
%     pix = pi.*(l-l_in-l_t); 
%     w = (sin(pix)./sin(pix./M).*exp(1i*(M-1)*pix/M))/M;    % equation 12 & 13
%     w(pix==0) = 1;
% end

function w = Sampling_Function_t(M,l,l_in,l_t)
    t = l-l_in-l_t;
    w = raised_cosine(t, 1, 0.5);
end