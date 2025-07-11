function [phi_trunc, phi_t, phi_t_aod, phi_t_aoa] = MIMO_Gen_measurement_matrix(x_p, x_kp, x_lp, Nr, Nt, M_T, N_T, M_aoa, N_aod, aod_bar, aoa_bar, delta_aod, delta_aoa)
% Massive MIMO �����������ɺ���
% ģ��ODDM�е�Gen_measurement_matrix
% ����:
%   x_p: ��Ƶ����
%   x_kp, x_lp: ��Ƶλ�ò���
%   Nr, Nt: ���պͷ���������
%   M_T, N_T: �ضϾ���ά��
%   M_aoa, N_aod: �Ƕ�����ά��
%   aod_bar, aoa_bar: �Ƕ�����
%   delta_aod, delta_aoa: �Ƕ�ƫ��
% ���:
%   phi_trunc: �ضϲ�������
%   phi_t: ������������
%   phi_t_aod: ����ǵ�������
%   phi_t_aoa: ���սǵ�������

phi_t = zeros(M_T * N_T, M_aoa * N_aod);
phi_t_aod = zeros(M_T * N_T, M_aoa * N_aod);
phi_t_aoa = zeros(M_T * N_T, M_aoa * N_aod);

for nn = 0:(N_aod-1)
    for mm = 1:M_aoa
        for nnn = 0:(N_T-1)
            for mmm = 1:M_T
                grid_idx = nn * M_aoa + mm;
                matrix_idx = nnn * M_T + mmm;
                
                % ������������
                phi_t(matrix_idx, grid_idx) = x_p * ...
                    MIMO_Array_Response_Tx(Nt, nnn, x_kp-1, aod_bar(grid_idx)) * ...
                    MIMO_Array_Response_Rx(Nr, mmm-1, x_lp-1, aoa_bar(grid_idx));
                
                % ����ǵ�������
                phi_t_aod(matrix_idx, grid_idx) = x_p * ...
                    MIMO_wt_derivation_aod(Nt, nnn, x_kp-1, aod_bar(grid_idx)) * ...
                    MIMO_Array_Response_Rx(Nr, mmm-1, x_lp-1, aoa_bar(grid_idx));
                
                % ���սǵ�������
                phi_t_aoa(matrix_idx, grid_idx) = x_p * ...
                    MIMO_Array_Response_Tx(Nt, nnn, x_kp-1, aod_bar(grid_idx)) * ...
                    MIMO_wt_derivation_aoa(Nr, mmm-1, x_lp-1, aoa_bar(grid_idx));
            end
        end
    end
end

phi_trunc = phi_t + phi_t_aod * diag(delta_aod) + phi_t_aoa * diag(delta_aoa);

end