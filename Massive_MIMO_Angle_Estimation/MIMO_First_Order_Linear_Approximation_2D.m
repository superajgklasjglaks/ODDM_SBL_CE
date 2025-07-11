function [aod_bar, aoa_bar] = MIMO_First_Order_Linear_Approximation_2D(N_aod, M_aoa, aod_max, r_aod, r_aoa)
% Massive MIMO 2Dһ�����Խ��ƺ���
% ģ��ODDM�е�First_Order_Linear_Approximation_2D
% ����:
%   N_aod: �����������
%   M_aoa: ���ս�������
%   aod_max: ��������ֵ
%   r_aod: ����Ƿֱ���
%   r_aoa: ���սǷֱ���
% ���:
%   aod_bar: �������������
%   aoa_bar: ���ս���������

% ���ɷ��������
aod_bar = linspace(-aod_max, aod_max, N_aod);
% ���ɽ��ս����� (������սǷ�Χ�뷢�����ͬ)
aoa_bar = linspace(-aod_max, aod_max, M_aoa);

end