function [h_hat, aod_hat, aoa_hat] = MIMO_hierarchical_SBL_refined_1D(x_p, x_kp, x_lp, Nr, Nt, N_T, M_T, y_T, refined_aod, refined_aoa, on_flag)
% Massive MIMO �ֲ�1D SBL ϸ������
% ģ��ODDM�е�hierarchical_SBL_refined_1D
% ����:
%   x_p: ��Ƶ����
%   x_kp, x_lp: ��Ƶλ�ò���
%   Nr, Nt: ���պͷ���������
%   N_T, M_T: �ضϾ���ά��
%   y_T: �����ź�����
%   refined_aod, refined_aoa: ϸ���ĽǶ�����
%   on_flag: �����־
% ���:
%   h_hat: ���Ƶ��ŵ�ϵ��
%   aod_hat, aoa_hat: ���ƵĽǶ�

virtual_size = length(refined_aod);
rho = 1e-2;
c = 1e-4;
d = 1e-4;

%% ����ϸ������Ĳ�������
delta_aod = zeros(virtual_size, 1);
delta_aoa = zeros(virtual_size, 1);

% ֱ��ʹ��ϸ��������Ϊ��������
aod_bar = refined_aod;
aoa_bar = refined_aoa;

% ���ɲ�������
phi_t = zeros(M_T * N_T, virtual_size);
for i = 1:virtual_size
    for nnn = 0:(N_T-1)
        for mmm = 1:M_T
            matrix_idx = nnn * M_T + mmm;
            phi_t(matrix_idx, i) = x_p * ...
                MIMO_Array_Response_Tx(Nt, nnn, x_kp-1, aod_bar(i)) * ...
                MIMO_Array_Response_Rx(Nr, mmm-1, x_lp-1, aoa_bar(i));
        end
    end
end

% �򻯵�SBL�㷨
ksi = 1e-3;
Tmax = 200;  % ���ٵ�������

sigma2_bar = sum(abs(y_T).^2) / (100 * N_T * M_T);
alpha = abs(phi_t.' * y_T);
beta0 = 1 / sigma2_bar;

% ��������
for t = 1:Tmax
    alpha_prev = alpha;
    delta_a = diag(alpha);
    
    % ����mu��Sigma
    C = (1 / beta0) * eye(M_T * N_T) + phi_t * delta_a * phi_t';
    C_inv = inv(C);
    sigma_hbar = delta_a - delta_a * phi_t' * C_inv * phi_t * delta_a;
    mu_hbar = beta0 * sigma_hbar * phi_t' * y_T;
    
    % ����gamma1
    gamma1 = 1 - real(diag(sigma_hbar)) ./ (alpha + 1e-16);
    
    % ����alpha
    alpha_temp = abs(mu_hbar).^2 + real(diag(sigma_hbar));
    alpha = -0.5/rho + sqrt(0.25 / rho^2 + alpha_temp / rho);
    
    % ����beta0
    resid = y_T - phi_t * mu_hbar;
    res2 = norm(resid, 'fro')^2;
    A_beta0 = res2 + 1 / beta0 * sum(gamma1);
    beta0 = (c - 1 + M_T * N_T) / (d + A_beta0);
    
    % ��ֹ����
    tolerance = norm(alpha - alpha_prev) / norm(alpha_prev);
    if tolerance <= ksi
        break;
    end
end

%% ���
h_hat = mu_hbar;
aod_hat = aod_bar;
aoa_hat = aoa_bar;

end