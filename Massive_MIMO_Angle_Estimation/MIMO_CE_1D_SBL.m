function [h_hat, aod_hat, aoa_hat, virtual_size, phi_trunc, delta_a] = MIMO_CE_1D_SBL(x_p, x_kp, x_lp, Nr, Nt, N_T, M_T, y_T, r_aod, r_aoa, aod_max, aoa_max, on_flag)
% Massive MIMO 1D SBL �Ƕȹ��ƺ���
% ģ��ODDM�е�CE_1D_SBL
% ����:
%   x_p: ��Ƶ����
%   x_kp, x_lp: ��Ƶλ�ò���
%   Nr, Nt: ���պͷ���������
%   N_T, M_T: �ضϾ���ά��
%   y_T: �����ź�����
%   r_aod, r_aoa: �Ƕȷֱ���
%   aod_max, aoa_max: �Ƕ����ֵ
%   on_flag: �����־
% ���:
%   h_hat: ���Ƶ��ŵ�ϵ��
%   aod_hat, aoa_hat: ���ƵĽǶ�
%   virtual_size: ���������С
%   phi_trunc: �ضϲ�������
%   delta_a: �ԽǾ���

N_aod = ceil(2 * aod_max / r_aod);   % ����������������
M_aoa = ceil(2 * aoa_max / r_aoa);   % ���ս������������
virtual_size = N_aod * M_aoa;
rho = 1e-2;
c = 1e-4;
d = 1e-4;

%% ��ʼ������Ƕȹ���
delta_aod = zeros(virtual_size, 1);
delta_aoa = zeros(virtual_size, 1);
[aod_bar, aoa_bar] = MIMO_First_Order_Linear_Approximation(N_aod, M_aoa, aod_max, r_aod, r_aoa);
[phi_trunc, phi_t, phi_t_aod, phi_t_aoa] = MIMO_Gen_measurement_matrix(x_p, x_kp, x_lp, Nr, Nt, M_T, N_T, M_aoa, N_aod, aod_bar, aoa_bar, delta_aod, delta_aoa);

ksi = 1e-3;
Tmax = 5e2;

sigma2_bar = sum(abs(y_T).^2) / (100 * N_T * M_T);
P_hat = floor(M_T * N_T / log(virtual_size));
alpha = abs(phi_trunc.' * y_T);
beta0 = 1 / sigma2_bar;

% ��������
for t = 1:Tmax
    phi_trunc = phi_t + phi_t_aod * diag(delta_aod) + phi_t_aoa * diag(delta_aoa);
    alpha_prev = alpha;
    delta_a = diag(alpha);
    
    % �㷨����3: ����mu��Sigma
    C = (1 / beta0) * eye(M_T * N_T) + phi_trunc * delta_a * phi_trunc';
    C_inv = inv(C);
    sigma_hbar = delta_a - delta_a * phi_trunc' * C_inv * phi_trunc * delta_a;  % ����41
    mu_hbar = beta0 * sigma_hbar * phi_trunc' * y_T;    % ����40
    
    % �㷨����4: ����gamma1���ڸ���beta0
    gamma1 = 1 - real(diag(sigma_hbar)) ./ (alpha + 1e-16);
    
    % ����alpha
    alpha_temp = abs(mu_hbar).^2 + real(diag(sigma_hbar));
    alpha = -0.5/rho + sqrt(0.25 / rho^2 + alpha_temp / rho);  % ����49
    
    % ����beta0
    resid = y_T - phi_trunc * mu_hbar;
    res2 = norm(resid, 'fro')^2;
    A_beta0 = res2 + 1 / beta0 * sum(gamma1); % ����51
    beta0 = (c - 1 + M_T * N_T) / (d + A_beta0);
    
    % ��ֹ����
    tolerance = norm(alpha - alpha_prev) / norm(alpha_prev);
    if tolerance <= ksi
        break;
    end
    
    [~, idx] = sort(alpha, 'descend');
    delta_aod_prev = delta_aod;
    delta_aoa_prev = delta_aoa;
    mmph = (mu_hbar * mu_hbar' + sigma_hbar);
    
    if (tolerance < 1000 * ksi)
        %% ����delta_aod (�����ƫ��)
        idx = idx(1:P_hat);
        
        pHpv = phi_t_aod' * phi_t_aod;
        tptlv = (phi_t + phi_t_aoa * diag(delta_aoa_prev));
        A_aod = real(pHpv .* (conj(mu_hbar) * mu_hbar.' + sigma_hbar.')); % ����54
        b_aod1 = diag(mu_hbar) * phi_t_aod.' * conj(y_T);
        b_aod2 = diag(mmph * tptlv' * phi_t_aod);
        b_aod = real(b_aod1 - b_aod2);    % ����55
        
        A_T_aod = A_aod(idx, idx);
        b_T_aod = b_aod(idx);
        
        temp = delta_aod;
        delta_T_aod = A_T_aod \ b_T_aod;
        delta_aod_upgrate = zeros(virtual_size, 1);
        delta_aod_upgrate(idx) = temp(idx);
        flag = ((max(svd(A_T_aod)) / min(svd(A_T_aod))) < 1e6);
        
        if (flag == 0)
            for pp = 1:P_hat
                temp_delta_aod = delta_aod(idx);
                temp_delta_aod(pp) = 0;
                if A_T_aod(pp, pp) == 0
                    delta_aod_upgrate(idx(pp)) = 0;
                    continue;
                else
                    delta_aod_upgrate(idx(pp)) = (b_T_aod(pp) - A_T_aod(pp, :) * temp_delta_aod) / A_T_aod(pp, pp);   % ����57
                end
            end
        else
            delta_aod_upgrate = zeros(virtual_size, 1);
            delta_aod_upgrate(idx) = delta_T_aod;
        end
        
        delta_aod_upgrate(delta_aod_upgrate < (-r_aod)/2) = (-r_aod)/2;
        delta_aod_upgrate(delta_aod_upgrate > (r_aod)/2) = (r_aod)/2;
        if (on_flag == 1)
            delta_aod = delta_aod_upgrate;
        end
        
        %% ����delta_aoa (���ս�ƫ��)
        pHpt = phi_t_aoa' * phi_t_aoa;
        tptlt = (phi_t + phi_t_aod * diag(delta_aod_prev));
        A_aoa = real(pHpt .* (conj(mu_hbar) * mu_hbar.' + sigma_hbar.')); % ����54
        b_aoa1 = diag(mu_hbar) * phi_t_aoa.' * conj(y_T);
        b_aoa2 = diag((mmph * tptlt') * phi_t_aoa);
        b_aoa = real(b_aoa1 - b_aoa2);    % ����55
        
        A_T_aoa = A_aoa(idx, idx);
        b_T_aoa = b_aoa(idx);
        
        temp = delta_aoa;
        delta_T_aoa = A_T_aoa \ b_T_aoa;
        delta_aoa_upgrate = zeros(virtual_size, 1);
        delta_aoa_upgrate(idx) = temp(idx);
        flag = ((max(svd(A_T_aoa)) / min(svd(A_T_aoa))) < 1e6);
        
        if (flag == 0)
            for pp = 1:P_hat
                temp_delta_aoa = delta_aoa(idx);
                temp_delta_aoa(pp) = 0;
                if A_T_aoa(pp, pp) == 0
                    delta_aoa_upgrate(idx(pp)) = 0;
                    continue;
                else
                    delta_aoa_upgrate(idx(pp)) = (b_T_aoa(pp) - A_T_aoa(pp, :) * temp_delta_aoa) / A_T_aoa(pp, pp);   % ����57
                end
            end
        else
            delta_aoa_upgrate = zeros(virtual_size, 1);
            delta_aoa_upgrate(idx) = delta_T_aoa;
        end
        
        delta_aoa_upgrate(delta_aoa_upgrate < (-r_aoa)/2) = (-r_aoa)/2;
        delta_aoa_upgrate(delta_aoa_upgrate > (r_aoa)/2) = (r_aoa)/2;
        if (on_flag == 1)
            delta_aoa = delta_aoa_upgrate;
        end
    end
end

%% ���
h_hat = mu_hbar;
aod_hat = aod_bar + delta_aod;
aoa_hat = aoa_bar + delta_aoa;

end