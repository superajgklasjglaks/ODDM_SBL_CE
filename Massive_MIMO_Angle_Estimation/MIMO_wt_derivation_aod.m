function w = MIMO_wt_derivation_aod(Nt, k, k_in, aod)
% Massive MIMO ����ǵ�������
% ģ��ODDM�е�wv_derivation
% ����:
%   Nt: ����������
%   k: ��ǰ��������
%   k_in: �ο���������
%   aod: ���䷽��� (����)
% ���:
%   w: �Ƕȵ���ֵ

% ���߼�� (�Բ���Ϊ��λ)
d_lambda = 0.5;

% ����Ƕȵ���
antenna_index = k - k_in;
w = 1j * 2 * pi * d_lambda * antenna_index * cos(aod) * ...
    exp(1j * 2 * pi * d_lambda * antenna_index * sin(aod)) / sqrt(Nt);

end