function [h_hat, omp_index, aod_hat, aoa_hat] = MIMO_OMP(x_p, x_kp, x_lp, Nr, Nt, N_T, M_T, y_T, r_aod, r_aoa, aod_max, aoa_max)
% Massive MIMO OMP �Ƕȹ��ƺ���
% ģ��ODDM�е�OMP
% ����:
%   x_p: ��Ƶ����
%   x_kp, x_lp: ��Ƶλ�ò���
%   Nr, Nt: ���պͷ���������
%   N_T, M_T: �ضϾ���ά��
%   y_T: �����ź�����
%   r_aod, r_aoa: �Ƕȷֱ���
%   aod_max, aoa_max: �Ƕ����ֵ
% ���:
%   h_hat: ���Ƶ��ŵ�ϵ��
%   omp_index: OMPѡ�������
%   aod_hat, aoa_hat: ���ƵĽǶ�

N_aod = ceil(2 * aod_max / r_aod);
M_aoa = ceil(2 * aoa_max / r_aoa);
virtual_size = N_aod * M_aoa;

% ���ɽǶ�����
[aod_bar, aoa_bar] = MIMO_First_Order_Linear_Approximation(N_aod, M_aoa, aod_max, r_aod, r_aoa);

% �����ֵ����
A = zeros(M_T * N_T, virtual_size);
for i = 1:virtual_size
    for nnn = 0:(N_T-1)
        for mmm = 1:M_T
            matrix_idx = nnn * M_T + mmm;
            A(matrix_idx, i) = x_p * ...
                MIMO_Array_Response_Tx(Nt, nnn, x_kp-1, aod_bar(i)) * ...
                MIMO_Array_Response_Rx(Nr, mmm-1, x_lp-1, aoa_bar(i));
        end
    end
end

% OMP�㷨����
max_iter = min(10, floor(M_T * N_T / 4));  % ����������
threshold = 1e-6;  % �в���ֵ

% ��ʼ��
residual = y_T;
selected_indices = [];
selected_atoms = [];

% OMP��ѭ��
for iter = 1:max_iter
    % �����ڻ�
    correlations = abs(A' * residual);
    
    % �ҵ��������Ե�ԭ��
    [~, max_idx] = max(correlations);
    
    % ��ӵ�ѡ�񼯺�
    selected_indices = [selected_indices, max_idx];
    selected_atoms = [selected_atoms, A(:, max_idx)];
    
    % ��С�������
    if size(selected_atoms, 2) == 1
        coeffs = selected_atoms \ y_T;
    else
        coeffs = pinv(selected_atoms) * y_T;
    end
    
    % ���²в�
    residual = y_T - selected_atoms * coeffs;
    
    % �����ֹ����
    if norm(residual) < threshold
        break;
    end
end

% ����ϡ���
h_hat = zeros(virtual_size, 1);
if ~isempty(selected_indices)
    h_hat(selected_indices) = coeffs;
end

% ���
omp_index = selected_indices;
aod_hat = aod_bar;
aoa_hat = aoa_bar;

end