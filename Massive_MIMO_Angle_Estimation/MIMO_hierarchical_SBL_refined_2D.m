function [H_opt, aod_opt, aoa_opt] = MIMO_hierarchical_SBL_refined_2D(x_p, x_kp, x_lp_factor, x_lp, Nr, Nt, N_T, M_T, y_trunc, aod_refined, aoa_refined, r_aod_fine, r_aoa_fine, aod_max, aoa_max, on_flag)
% Massive MIMO �ֲ�2D SBL ϸ������
% ģ��ODDM�е�hierarchical_SBL_refined_2D
% ����:
%   x_p: ��Ƶ����
%   x_kp, x_lp_factor, x_lp: ��Ƶλ�ò���
%   Nr, Nt: ���պͷ���������
%   N_T, M_T: �ضϾ���ά��
%   y_trunc: �ضϽ����źž���
%   aod_refined, aoa_refined: ϸ���ĽǶ�����
%   r_aod_fine, r_aoa_fine: ϸ����ֱ���
%   aod_max, aoa_max: �Ƕ����ֵ
%   on_flag: �����־
% ���:
%   H_opt: ���Ƶ�2D�ŵ�����
%   aod_opt, aoa_opt: ���ƵĽǶ�����

% ��y_trunc��������Ϊ����
y_T = reshape(y_trunc.', N_T * M_T, 1);

% ���÷ֲ�1D SBL��Ϊ����
[h_hat_temp, aod_hat_temp, aoa_hat_temp] = MIMO_hierarchical_SBL_refined_1D(x_p, x_kp, x_lp, Nr, Nt, N_T, M_T, y_T, aod_refined, aoa_refined, on_flag);

% ��������ά��
N_aod_refined = length(aod_refined);
M_aoa_refined = length(aoa_refined);

% �������������Ϊ2D��ʽ
if N_aod_refined * M_aoa_refined == length(h_hat_temp)
    H_opt = reshape(h_hat_temp, M_aoa_refined, N_aod_refined);
else
    % ���ά�Ȳ�ƥ�䣬����һ�����ʴ�С�ľ���
    H_opt = zeros(M_aoa_refined, N_aod_refined);
    % ������ϵ�����ں��ʵ�λ��
    [~, max_idx] = max(abs(h_hat_temp));
    if ~isempty(max_idx)
        [row_idx, col_idx] = ind2sub([M_aoa_refined, N_aod_refined], max_idx);
        H_opt(row_idx, col_idx) = h_hat_temp(max_idx);
    end
end

aod_opt = aod_hat_temp;
aoa_opt = aoa_hat_temp;

end