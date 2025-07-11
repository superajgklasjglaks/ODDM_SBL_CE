function [h_hat, aod_hat, aoa_hat] = MIMO_traditional_beamforming(x_p, y_trunc, aod_range, aoa_range, sigma_2_p)
% Massive MIMO ��ͳ�������νǶȹ��ƺ���
% ģ��ODDM�е�traditional_impulse
% ����:
%   x_p: ��Ƶ����
%   y_trunc: �ضϽ����źž���
%   aod_range, aoa_range: �Ƕ�������Χ
%   sigma_2_p: ��������
% ���:
%   h_hat: ���Ƶ��ŵ�ϵ��
%   aod_hat, aoa_hat: ���ƵĽǶ�

[M_T, N_T] = size(y_trunc);

% �򵥵���������
aod_grid = linspace(-pi/3, pi/3, 2*aod_range+1);
aoa_grid = linspace(-pi/3, pi/3, 2*aoa_range+1);

virtual_size = length(aod_grid) * length(aoa_grid);
h_hat = zeros(virtual_size, 1);
aod_hat = zeros(virtual_size, 1);
aoa_hat = zeros(virtual_size, 1);

% ���߼��
d_lambda = 0.5;

index = 1;
for i = 1:length(aod_grid)
    for j = 1:length(aoa_grid)
        aod = aod_grid(i);
        aoa = aoa_grid(j);
        
        % ������������
        at = exp(1j * 2 * pi * d_lambda * (0:N_T-1)' * sin(aod)) / sqrt(N_T);
        ar = exp(1j * 2 * pi * d_lambda * (0:M_T-1)' * sin(aoa)) / sqrt(M_T);
        
        % ��������
        beamformed_signal = ar' * y_trunc * at;
        
        % �����ŵ�ϵ��
        h_hat(index) = beamformed_signal / x_p;
        aod_hat(index) = aod;
        aoa_hat(index) = aoa;
        
        index = index + 1;
    end
end

% �����������
noise_threshold = sqrt(sigma_2_p) * 3;
h_hat(abs(h_hat) < noise_threshold) = 0;

end