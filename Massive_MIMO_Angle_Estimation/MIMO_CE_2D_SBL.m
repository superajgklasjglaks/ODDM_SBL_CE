function [H_opt, aod_opt, aoa_opt, virtual_size] = MIMO_CE_2D_SBL(x_p, x_kp, x_lp_factor, x_lp, Nr, Nt, N_T, M_T, y_trunc, r_aod, r_aoa, aod_max, aoa_max, on_flag)
% Massive MIMO 2D SBL �Ƕȹ��ƺ���
% �򻯰汾��ģ��ODDM�е�CE_2D_SBL
% ����:
%   x_p: ��Ƶ����
%   x_kp, x_lp_factor, x_lp: ��Ƶλ�ò���
%   Nr, Nt: ���պͷ���������
%   N_T, M_T: �ضϾ���ά��
%   y_trunc: �ضϽ����źž���
%   r_aod, r_aoa: �Ƕȷֱ���
%   aod_max, aoa_max: �Ƕ����ֵ
%   on_flag: �����־
% ���:
%   H_opt: ���Ƶ�2D�ŵ�����
%   aod_opt, aoa_opt: ���ƵĽǶ�����
%   virtual_size: ���������С

N_aod = ceil(2 * aod_max / r_aod);
M_aoa = ceil(2 * aoa_max / r_aoa);
virtual_size = N_aod * M_aoa;

% ���ɽǶ�����
[aod_bar, aoa_bar] = MIMO_First_Order_Linear_Approximation_2D(N_aod, M_aoa, aod_max, r_aod, r_aoa);

% �򻯵�2D SBLʵ��
% ��y_trunc��������Ϊ����
y_T = reshape(y_trunc.', N_T * M_T, 1);

% ����1D SBL��Ϊ����
[h_hat_temp, aod_hat_temp, aoa_hat_temp, ~, ~, ~] = MIMO_CE_1D_SBL(x_p, x_kp, x_lp, Nr, Nt, N_T, M_T, y_T, r_aod, r_aoa, aod_max, aoa_max, on_flag);

% �������������Ϊ2D��ʽ
H_opt = reshape(h_hat_temp, M_aoa, N_aod);
aod_opt = aod_hat_temp;
aoa_opt = aoa_hat_temp;

end