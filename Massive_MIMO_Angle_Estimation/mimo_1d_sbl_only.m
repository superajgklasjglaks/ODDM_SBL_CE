% Massive MIMO 1D SBL ����ǹ����㷨
% �򻯰汾��ֻ����1D SBL�㷨
% ����: AI Assistant
% ����: 2024

clc; clear; close all;

%% ϵͳ��������
% Massive MIMOϵͳ����
Nt = 32;                    % ����������
Nr = 32;                    % ����������
P = 3;                      % ·����
L_pilot = 16;               % ��Ƶ����

% �ǶȲ���
aod_range = [-pi/3, pi/3];  % ����Ƿ�Χ (����)
aoa_range = [-pi/3, pi/3];  % ���սǷ�Χ (����)

% �㷨����
resolution_1D = 0.1;       % 1D����ֱ���
virtual_factor = 4;         % ������������

% �������
num_frames = 100;           % ����֡��
SNR_dB = 10;                % ����� (dB)

%% ������ʵ�Ƕ�
% ������ɷ���Ǻͽ��ս�
aod_true = aod_range(1) + (aod_range(2) - aod_range(1)) * rand(P, 1);
aoa_true = aoa_range(1) + (aoa_range(2) - aoa_range(1)) * rand(P, 1);

% ������ʵ�ŵ�ϵ��
h_true = (randn(P, 1) + 1j * randn(P, 1)) / sqrt(2);

fprintf('��ʵ����� (��): ');
fprintf('%.2f ', aod_true * 180 / pi);
fprintf('\n');
fprintf('��ʵ���ս� (��): ');
fprintf('%.2f ', aoa_true * 180 / pi);
fprintf('\n\n');

%% 1D SBL�㷨����
% ���ɽǶ�����
aod_grid = aod_range(1):resolution_1D:aod_range(2);
aoa_grid = aoa_range(1):resolution_1D:aoa_range(2);
G_aod = length(aod_grid);
G_aoa = length(aoa_grid);

fprintf('������������: %d\n', G_aod);
fprintf('���ս��������: %d\n', G_aoa);

%% ��ʼ������洢
NMSE_aod_1D = zeros(num_frames, 1);
NMSE_aoa_1D = zeros(num_frames, 1);
NMSE_h_1D = zeros(num_frames, 1);

%% Monte Carlo����
fprintf('\n��ʼ1D SBL�㷨����...\n');
tic;

for frame = 1:num_frames
    if mod(frame, 20) == 0
        fprintf('���֡��: %d/%d\n', frame, num_frames);
    end
    
    %% ���ɵ�Ƶ�źźͽ����ź�
    % ��Ƶ���� (�����Ƶ)
    X_pilot = (randn(Nt, L_pilot) + 1j * randn(Nt, L_pilot)) / sqrt(2);
    
    % �����ŵ�����
    H = zeros(Nr, Nt);
    for p = 1:P
        % ���䵼������
        a_tx = exp(1j * (0:Nt-1)' * pi * sin(aod_true(p)));
        % ���յ�������
        a_rx = exp(1j * (0:Nr-1)' * pi * sin(aoa_true(p)));
        % �ŵ�����
        H = H + h_true(p) * a_rx * a_tx';
    end
    
    % �����ź�
    Y_pilot = H * X_pilot;
    
    % �������
    noise_power = 10^(-SNR_dB/10);
    noise = sqrt(noise_power/2) * (randn(Nr, L_pilot) + 1j * randn(Nr, L_pilot));
    Y_pilot_noisy = Y_pilot + noise;
    
    % �����������ź�
    y = Y_pilot_noisy(:);
    
    %% 1D SBL�㷨
    try
        [h_hat_1D, aod_hat_1D, aoa_hat_1D] = MIMO_CE_1D_SBL(...
            y, X_pilot, Nt, Nr, L_pilot, P, ...
            aod_grid, aoa_grid, virtual_factor);
        
        % ����NMSE
        if ~isempty(aod_hat_1D) && ~isempty(aoa_hat_1D)
            % �Ƕ�ƥ�� (�ҵ���ӽ��Ĺ���ֵ)
            aod_error = min(abs(aod_hat_1D - aod_true'));
            aoa_error = min(abs(aoa_hat_1D - aoa_true'));
            
            NMSE_aod_1D(frame) = mean(aod_error.^2) / mean(aod_true.^2);
            NMSE_aoa_1D(frame) = mean(aoa_error.^2) / mean(aoa_true.^2);
            
            if ~isempty(h_hat_1D)
                h_error = min(abs(h_hat_1D - h_true'));
                NMSE_h_1D(frame) = mean(h_error.^2) / mean(abs(h_true).^2);
            else
                NMSE_h_1D(frame) = 1;
            end
        else
            NMSE_aod_1D(frame) = 1;
            NMSE_aoa_1D(frame) = 1;
            NMSE_h_1D(frame) = 1;
        end
        
    catch ME
        fprintf('��%d֡1D SBL�㷨����: %s\n', frame, ME.message);
        NMSE_aod_1D(frame) = 1;
        NMSE_aoa_1D(frame) = 1;
        NMSE_h_1D(frame) = 1;
    end
end

simulation_time = toc;
fprintf('������ɣ��ܺ�ʱ: %.2f��\n\n', simulation_time);

%% ���ͳ��
fprintf('=== 1D SBL�㷨����ͳ�� ===\n');
fprintf('�����NMSE (dB): %.2f\n', 10*log10(mean(NMSE_aod_1D)));
fprintf('���ս�NMSE (dB): %.2f\n', 10*log10(mean(NMSE_aoa_1D)));
fprintf('�ŵ�ϵ��NMSE (dB): %.2f\n', 10*log10(mean(NMSE_h_1D)));
fprintf('ƽ��ÿ֡��ʱ: %.4f��\n', simulation_time/num_frames);

%% ���ƽ��
figure('Position', [100, 100, 1200, 400]);

% �����NMSE
subplot(1, 3, 1);
plot(1:num_frames, 10*log10(NMSE_aod_1D), 'b-', 'LineWidth', 1.5);
hold on;
yline(10*log10(mean(NMSE_aod_1D)), 'r--', 'LineWidth', 2);
xlabel('����֡��');
ylabel('�����NMSE (dB)');
title('1D SBL - ����ǹ�������');
grid on;
legend('˲ʱNMSE', sprintf('ƽ��NMSE: %.2f dB', 10*log10(mean(NMSE_aod_1D))), 'Location', 'best');

% ���ս�NMSE
subplot(1, 3, 2);
plot(1:num_frames, 10*log10(NMSE_aoa_1D), 'g-', 'LineWidth', 1.5);
hold on;
yline(10*log10(mean(NMSE_aoa_1D)), 'r--', 'LineWidth', 2);
xlabel('����֡��');
ylabel('���ս�NMSE (dB)');
title('1D SBL - ���սǹ�������');
grid on;
legend('˲ʱNMSE', sprintf('ƽ��NMSE: %.2f dB', 10*log10(mean(NMSE_aoa_1D))), 'Location', 'best');

% �ŵ�ϵ��NMSE
subplot(1, 3, 3);
plot(1:num_frames, 10*log10(NMSE_h_1D), 'm-', 'LineWidth', 1.5);
hold on;
yline(10*log10(mean(NMSE_h_1D)), 'r--', 'LineWidth', 2);
xlabel('����֡��');
ylabel('�ŵ�ϵ��NMSE (dB)');
title('1D SBL - �ŵ�ϵ����������');
grid on;
legend('˲ʱNMSE', sprintf('ƽ��NMSE: %.2f dB', 10*log10(mean(NMSE_h_1D))), 'Location', 'best');

sgtitle('Massive MIMO 1D SBL ����ǹ�������', 'FontSize', 14, 'FontWeight', 'bold');

%% ��ʾ���չ��ƽ��
fprintf('\n=== ���һ֡�Ĺ��ƽ�� ===\n');
if exist('aod_hat_1D', 'var') && ~isempty(aod_hat_1D)
    fprintf('���Ʒ���� (��): ');
    fprintf('%.2f ', aod_hat_1D * 180 / pi);
    fprintf('\n');
end
if exist('aoa_hat_1D', 'var') && ~isempty(aoa_hat_1D)
    fprintf('���ƽ��ս� (��): ');
    fprintf('%.2f ', aoa_hat_1D * 180 / pi);
    fprintf('\n');
end
if exist('h_hat_1D', 'var') && ~isempty(h_hat_1D)
    fprintf('�����ŵ�ϵ������: ');
    fprintf('%.4f ', abs(h_hat_1D));
    fprintf('\n');
end

fprintf('\n1D SBL�㷨������ɣ�\n');