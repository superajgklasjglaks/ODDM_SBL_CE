close all
clear all
clc

%% 绘制算法性能随路径数量变化的图表
% 基于数值结果分析表格数据

%% 数据定义
% 路径数量
P_values = [3, 6, 9, 12, 15];

% NMSE性能数据 (dB)
NMSE_1D_SBL = [-23.16, -24.34, -22.45, -22.15, -21.57];
NMSE_2D_SBL = [-22.92, -21.81, -18.76, -16.78, -16.01];
NMSE_hierarchical_1D = [-23.27, -24.03, -21.79, -21.06, -20.95];
NMSE_hierarchical_2D = [-0.47, -0.52, -0.66, -0.47, -0.12];
NMSE_OMP = [-21.28, -20.71, -19.14, -18.05, -17.70];
NMSE_Traditional = [-15.29, -14.13, -12.53, -13.95, -14.86];

%% 颜色定义
colors = [0,107,182;     % 蓝色 - 1D SBL
          118,41,133;    % 紫色 - 2D SBL
          234,174,31;    % 黄色 - 分层1D SBL
          215,94,59;     % 橙色 - 分层2D SBL
          184,125,182;   % 淡紫色 - OMP
          71,90,40]/256; % 绿色 - Traditional

%% 创建图表
figure('Position', [100, 100, 1200, 800]);

%% 子图1: 主要算法NMSE性能对比
plot(P_values, NMSE_1D_SBL, '-s', 'Color', colors(1,:), 'LineWidth', 2.5, 'MarkerSize', 8, 'MarkerFaceColor', colors(1,:));
hold on;
plot(P_values, NMSE_2D_SBL, '-o', 'Color', colors(2,:), 'LineWidth', 2.5, 'MarkerSize', 8, 'MarkerFaceColor', colors(2,:));
plot(P_values, NMSE_hierarchical_1D, '-^', 'Color', colors(3,:), 'LineWidth', 2.5, 'MarkerSize', 8, 'MarkerFaceColor', colors(3,:));
plot(P_values, NMSE_OMP, '-d', 'Color', colors(5,:), 'LineWidth', 2.5, 'MarkerSize', 8, 'MarkerFaceColor', colors(5,:));
plot(P_values, NMSE_Traditional, '-p', 'Color', colors(6,:), 'LineWidth', 2.5, 'MarkerSize', 8, 'MarkerFaceColor', colors(6,:));
grid on;
xlabel('路径数量 P', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('NMSE (dB)', 'FontSize', 12, 'FontWeight', 'bold');
title('主要算法NMSE性能随路径数量变化', 'FontSize', 14, 'FontWeight', 'bold');
legend('1D SBL', '2D SBL', '分层1D SBL', 'OMP', 'Traditional', 'Location', 'best', 'FontSize', 10);
set(gca, 'FontSize', 11);
ylim([-26, -10]);
