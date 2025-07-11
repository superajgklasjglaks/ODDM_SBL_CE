function w = MIMO_Array_Response_Rx(Nr, l, l_in, aoa)
% Massive MIMO ����������Ӧ����
% ģ��ODDM�е�Sampling_Function_t
% ����:
%   Nr: ����������
%   l: ��ǰ��������
%   l_in: �ο���������
%   aoa: ���շ���� (����)
% ���:
%   w: ������Ӧֵ

% ���߼�� (�Բ���Ϊ��λ)
d_lambda = 0.5;

% ����������Ӧ
% ʹ�þ�����������ģ��
antenna_index = l - l_in;
w = exp(1j * 2 * pi * d_lambda * antenna_index * sin(aoa)) / sqrt(Nr);

end