function [aod_bar, aoa_bar] = MIMO_First_Order_Linear_Approximation(N_aod, M_aoa, aod_max, r_aod, r_aoa)
% Massive MIMO һ�����Խ��ƺ���
% ģ��ODDM�е�First_Order_Linear_Approximation
% ����:
%   N_aod: �����������
%   M_aoa: ���ս�������
%   aod_max: ��������ֵ
%   r_aod: ����Ƿֱ���
%   r_aoa: ���սǷֱ���
% ���:
%   aod_bar: ���������
%   aoa_bar: ���ս�����

virtual_size = N_aod * M_aoa;
aod_bar = zeros(virtual_size, 1);
aoa_bar = zeros(virtual_size, 1);

% ���ɷ��������
aod_grid = linspace(-aod_max, aod_max, N_aod);
% ���ɽ��ս�����
aoa_grid = linspace(-aod_max, aod_max, M_aoa);  % ������սǷ�Χ�뷢�����ͬ

% �������
index = 1;
for n = 1:N_aod
    for m = 1:M_aoa
        aod_bar(index) = aod_grid(n);
        aoa_bar(index) = aoa_grid(m);
        index = index + 1;
    end
end

end