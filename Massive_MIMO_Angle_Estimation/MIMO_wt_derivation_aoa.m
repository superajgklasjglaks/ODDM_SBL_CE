function w = MIMO_wt_derivation_aoa(Nr, l, l_in, aoa)
% Massive MIMO ���սǵ�������
% ģ��ODDM�е�wt_derivation
% ����:
%   Nr: ����������
%   l: ��ǰ��������
%   l_in: �ο���������
%   aoa: ���շ���� (����)
% ���:
%   w: �Ƕȵ���ֵ

% ���߼�� (�Բ���Ϊ��λ)
d_lambda = 0.5;

% ����Ƕȵ���
antenna_index = l - l_in;
w = 1j * 2 * pi * d_lambda * antenna_index * cos(aoa) * ...
    exp(1j * 2 * pi * d_lambda * antenna_index * sin(aoa)) / sqrt(Nr);

end