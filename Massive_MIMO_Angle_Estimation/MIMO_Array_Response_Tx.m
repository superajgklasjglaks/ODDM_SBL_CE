function w = MIMO_Array_Response_Tx(Nt, k, k_in, aod)
% Massive MIMO ����������Ӧ����
% ģ��ODDM�е�Sampling_Function_v
% ����:
%   Nt: ����������
%   k: ��ǰ��������
%   k_in: �ο���������
%   aod: ���䷽��� (����)
% ���:
%   w: ������Ӧֵ

% ���߼�� (�Բ���Ϊ��λ)
d_lambda = 0.5;

% ����������Ӧ
% ʹ�þ�����������ģ��
antenna_index = k - k_in;
w = exp(1j * 2 * pi * d_lambda * antenna_index * sin(aod)) / sqrt(Nt);

end