close all
clear all
clc

%% �����㷨������·�������仯��ͼ��
% ������ֵ��������������

%% ���ݶ���
% ·������
P_values = [3, 6, 9, 12, 15];

% NMSE�������� (dB)
NMSE_1D_SBL = [-23.16, -24.34, -22.45, -22.15, -21.57];
NMSE_2D_SBL = [-22.92, -21.81, -18.76, -16.78, -16.01];
NMSE_hierarchical_1D = [-23.27, -24.03, -21.79, -21.06, -20.95];
NMSE_hierarchical_2D = [-0.47, -0.52, -0.66, -0.47, -0.12];
NMSE_OMP = [-21.28, -20.71, -19.14, -18.05, -17.70];
NMSE_Traditional = [-15.29, -14.13, -12.53, -13.95, -14.86];

%% ��ɫ����
colors = [0,107,182;     % ��ɫ - 1D SBL
          118,41,133;    % ��ɫ - 2D SBL
          234,174,31;    % ��ɫ - �ֲ�1D SBL
          215,94,59;     % ��ɫ - �ֲ�2D SBL
          184,125,182;   % ����ɫ - OMP
          71,90,40]/256; % ��ɫ - Traditional

%% ����ͼ��
figure('Position', [100, 100, 1200, 800]);

%% ��ͼ1: ��Ҫ�㷨NMSE���ܶԱ�
plot(P_values, NMSE_1D_SBL, '-s', 'Color', colors(1,:), 'LineWidth', 2.5, 'MarkerSize', 8, 'MarkerFaceColor', colors(1,:));
hold on;
plot(P_values, NMSE_2D_SBL, '-o', 'Color', colors(2,:), 'LineWidth', 2.5, 'MarkerSize', 8, 'MarkerFaceColor', colors(2,:));
plot(P_values, NMSE_hierarchical_1D, '-^', 'Color', colors(3,:), 'LineWidth', 2.5, 'MarkerSize', 8, 'MarkerFaceColor', colors(3,:));
plot(P_values, NMSE_OMP, '-d', 'Color', colors(5,:), 'LineWidth', 2.5, 'MarkerSize', 8, 'MarkerFaceColor', colors(5,:));
plot(P_values, NMSE_Traditional, '-p', 'Color', colors(6,:), 'LineWidth', 2.5, 'MarkerSize', 8, 'MarkerFaceColor', colors(6,:));
grid on;
xlabel('·������ P', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('NMSE (dB)', 'FontSize', 12, 'FontWeight', 'bold');
title('��Ҫ�㷨NMSE������·�������仯', 'FontSize', 14, 'FontWeight', 'bold');
legend('1D SBL', '2D SBL', '�ֲ�1D SBL', 'OMP', 'Traditional', 'Location', 'best', 'FontSize', 10);
set(gca, 'FontSize', 11);
ylim([-26, -10]);
