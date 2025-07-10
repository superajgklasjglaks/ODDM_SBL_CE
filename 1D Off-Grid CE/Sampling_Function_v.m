function w = Sampling_Function_v(N,k,k_in,k_v)
    pix = pi.*(k-k_in-k_v); 
    w = (sin(pix)./sin(pix./N).*exp(-1i*(N-1).*pix./N))/N;    % equation 12 & 13
    w(pix==0) = 1;
end

%% Sampling_Function(32,1,1,0)
%% Sampling_Function(32,32,32,0)
%{
a = 0:0.02:19;
N = 20;
H1 = Sampling_Function(N,a,6,0);
H = abs(Sampling_Function(N,a,6,0));
figure;
plot(a,H);

%}