close all
clear all
clc

%% SNR
SNR_dB = 0:5:20;

NMSE_dB_trd = [-7.424, -13.552, -14.966, -15.177, -15.182];
NMSE_dB_1D = [-7.994, -18.376, -26.492, -31.167, -31.689];
NMSE_dB_4D = [-11.770, -19.684, -23.819, -25.863, -26.451];
NMSE_dB_omp = [-8.137, -16.784, -20.949, -22.334, -22.519];
NMSE_dB_4Dsomp = [-8.937, -17.684, -22.249, -24.334, -25.119];
NMSE_dB_CRLB = [-20.015, -25.015, -30.015, -35.015, -40.015];

colors = [0,107,182; %蓝1
          118,41,133; %紫2
          234,174,31; %黄3
          215,94,59; %橙4
          184,125,182; %粉紫5
          71,90,40; %绿6
          161,27,30]; %红7
colors = colors/256;

figure()

plot(SNR_dB,NMSE_dB_trd,'-s', 'Color', colors(1,:),'LineWidth',1.5,'MarkerSize',8);
hold on
plot(SNR_dB,NMSE_dB_1D,'-v', 'Color', colors(3,:),'LineWidth',1.5,'MarkerSize',8);
hold on
plot(SNR_dB,NMSE_dB_4D,'-*', 'Color', colors(2,:),'LineWidth',1.5,'MarkerSize',8);
hold on
plot(SNR_dB,NMSE_dB_omp,'-o', 'Color', colors(5,:),'LineWidth',1.5,'MarkerSize',8);
hold on
plot(SNR_dB,NMSE_dB_4Dsomp,'-^', 'Color', colors(6,:),'LineWidth',1.5,'MarkerSize',8);
% hold on
% plot(SNR_dB,NMSE_dB_CRLB,'-x', 'Color', colors(4,:),'LineWidth',1.5,'MarkerSize',8);
grid on

legend('Traditional impulse','1D SBL', '4D SBL','OMP','4D-SOMP');
xlabel('SNR(dB)')
ylabel('NMSE(dB)')

%% 分辨率
virtual = 0.2:0.2:1;

NMSE_dB_trd = [-15.382, -14.312, -14.102, -13.831, -13.959];
NMSE_dB_1D = [-32.705, -27.763, -20.856, -17.328, -11.202];
NMSE_dB_4D = [-28.911, -24.938, -19.034, -16.646, -10.332];
NMSE_dB_omp = [-25.816, -20.422, -15.934, -13.660, -11.521];
NMSE_dB_4Dsomp = [-27.937, -22.684, -17.249, -15.834, -12.119];
%NMSE_dB_CRLB = [0.214, 0.068, 0.021, 0.007, 0.002];

colors = [0,107,182; %蓝1
          118,41,133; %紫2
          234,174,31; %黄3
          215,94,59; %橙4
          184,125,182; %粉紫5
          71,90,40; %绿6
          161,27,30]; %红7
colors = colors/256;

figure()

plot(virtual,NMSE_dB_trd,'-s', 'Color', colors(1,:),'LineWidth',1.5,'MarkerSize',8);
hold on
plot(virtual,NMSE_dB_1D,'-v', 'Color', colors(3,:),'LineWidth',1.5,'MarkerSize',8);
hold on
plot(virtual,NMSE_dB_4D,'-*', 'Color', colors(2,:),'LineWidth',1.5,'MarkerSize',8);
hold on
plot(virtual,NMSE_dB_omp,'-o', 'Color', colors(5,:),'LineWidth',1.5,'MarkerSize',8);
hold on
plot(virtual,NMSE_dB_4Dsomp,'-^', 'Color', colors(6,:),'LineWidth',1.5,'MarkerSize',8);
grid on

legend('Traditional impulse','1D SBL', '4D SBL','OMP','4D-SOMP');
xlabel('Virtual Resolution')
ylabel('NMSE(dB)')

%% 速度
speed = 200:200:1000;

NMSE_dB_trd = [-15.382, -14.312, -14.302, -13.831, -13.159];
NMSE_dB_1D = [-26.963, -26.763, -26.163, -26.563, -25.863];
NMSE_dB_4D = [-24.138, -24.938, -24.238, -25.138, -24.638];
NMSE_dB_omp = [-22.122, -22.422, -22.922, -22.822, -22.322];
NMSE_dB_4Dsomp = [-23.684, -23.384, -23.184, -22.584, -22.184];

colors = [0,107,182; %蓝1
          118,41,133; %紫2
          234,174,31; %黄3
          215,94,59; %橙4
          184,125,182; %粉紫5
          71,90,40; %绿6
          161,27,30]; %红7
colors = colors/256;

figure()

plot(speed,NMSE_dB_trd,'-s', 'Color', colors(1,:),'LineWidth',1.5,'MarkerSize',8);
hold on
plot(speed,NMSE_dB_1D,'-v', 'Color', colors(3,:),'LineWidth',1.5,'MarkerSize',8);
hold on
plot(speed,NMSE_dB_4D,'-*', 'Color', colors(2,:),'LineWidth',1.5,'MarkerSize',8);
hold on
plot(speed,NMSE_dB_omp,'-o', 'Color', colors(5,:),'LineWidth',1.5,'MarkerSize',8);
hold on
plot(speed,NMSE_dB_4Dsomp,'-^', 'Color', colors(6,:),'LineWidth',1.5,'MarkerSize',8);
grid on

legend('Traditional impulse','1D SBL', '4D SBL','OMP','4D-SOMP');
xlabel('User Velocity(m/s)')
ylabel('NMSE(dB)')

%% 天线数
% 定义横轴数据（数值）
x_values = [0, 8, 16, 24, 32];  

% 定义横轴标签（字符串）
x_labels = {'1×1', '8×8', '16×16', '24×24', '32×32'};  

NMSE_dB_trd = [-13.159, -13.831, -14.302, -14.312, -14.382];
NMSE_dB_1D = [-25.963, -26.763, -27.163, -27.563, -27.663];
NMSE_dB_4D = [-21.138, -22.938, -23.738, -24.138, -24.338];
NMSE_dB_omp = [-18.622, -20.322, -21.522, -22.322, -22.822];
NMSE_dB_4Dsomp = [-19.684, -21.784, -22.684, -23.184, -23.684];

colors = [0,107,182; %蓝1
          118,41,133; %紫2
          234,174,31; %黄3
          215,94,59; %橙4
          184,125,182; %粉紫5
          71,90,40; %绿6
          161,27,30]; %红7
colors = colors/256;

figure()

plot(x_values,NMSE_dB_trd,'-s', 'Color', colors(1,:),'LineWidth',1.5,'MarkerSize',8);
hold on
plot(x_values,NMSE_dB_1D,'-v', 'Color', colors(3,:),'LineWidth',1.5,'MarkerSize',8);
hold on
plot(x_values,NMSE_dB_4D,'-*', 'Color', colors(2,:),'LineWidth',1.5,'MarkerSize',8);
hold on
plot(x_values,NMSE_dB_omp,'-o', 'Color', colors(5,:),'LineWidth',1.5,'MarkerSize',8);
hold on
plot(x_values,NMSE_dB_4Dsomp,'-^', 'Color', colors(6,:),'LineWidth',1.5,'MarkerSize',8);
grid on

% 设置横轴标签
xticks(x_values);          % 设置刻度位置
xticklabels(x_labels);     % 设置刻度标签
xlim([0, 32]);

% 设置图形样式
grid on;
legend('Traditional impulse','1D SBL', '4D SBL','OMP','4D-SOMP');
xlabel('Antenna Configuration (N_r×N_t)');
ylabel('NMSE (dB)');

%% 28GHz and 625kHz
SNR_dB = 0:5:20;

NMSE_dB_trd = [-7.019, -12.424, -13.477, -13.606, -13.635];
NMSE_dB_1D = [-7.831, -17.887, -26.356, -28.925, -31.449];
NMSE_dB_4D = [-9.986, -17.767, -23.457, -25.470, -27.397];
NMSE_dB_omp = [-8.355, -17.330, -21.900, -23.332, -23.565];
NMSE_dB_4Dsomp = [-8.937, -17.684, -22.249, -24.334, -25.119];
NMSE_dB_CRLB = [-20.015, -25.015, -30.015, -35.015, -40.015];

colors = [0,107,182; %蓝1
          118,41,133; %紫2
          234,174,31; %黄3
          215,94,59; %橙4
          184,125,182; %粉紫5
          71,90,40; %绿6
          161,27,30]; %红7
colors = colors/256;

figure()

plot(SNR_dB,NMSE_dB_trd,'-s', 'Color', colors(1,:),'LineWidth',1.5,'MarkerSize',8);
hold on
plot(SNR_dB,NMSE_dB_1D,'-v', 'Color', colors(3,:),'LineWidth',1.5,'MarkerSize',8);
hold on
plot(SNR_dB,NMSE_dB_4D,'-*', 'Color', colors(2,:),'LineWidth',1.5,'MarkerSize',8);
hold on
plot(SNR_dB,NMSE_dB_omp,'-o', 'Color', colors(5,:),'LineWidth',1.5,'MarkerSize',8);
hold on
plot(SNR_dB,NMSE_dB_4Dsomp,'-^', 'Color', colors(6,:),'LineWidth',1.5,'MarkerSize',8);
hold on
% plot(SNR_dB,NMSE_dB_CRLB,'-x', 'Color', colors(4,:),'LineWidth',1.5,'MarkerSize',8);
% grid on

legend('Traditional impulse','1D SBL', '4D SBL','OMP','4D-SOMP');
xlabel('SNR(dB)')
ylabel('NMSE(dB)')

%% 导频数量
SNR_dB = 1:2:9;

NMSE_dB_trd = [-15.146, -15.270, -16.363, -16.966, -17.389];
NMSE_dB_1D = [-23.667, -25.492, -27.181, -28.331, -28.458];
NMSE_dB_4D = [-21.653, -22.920, -23.819, -24.307, -24.513];
NMSE_dB_omp = [-18.220, -19.902, -21.945, -22.959, -23.249];
NMSE_dB_4Dsomp = [-20.220, -22.102, -23.345, -23.959, -24.249];

colors = [0,107,182; %蓝1
          118,41,133; %紫2
          234,174,31; %黄3
          215,94,59; %橙4
          184,125,182; %粉紫5
          71,90,40; %绿6
          161,27,30]; %红7
colors = colors/256;

figure()

plot(SNR_dB,NMSE_dB_trd,'-s', 'Color', colors(1,:),'LineWidth',1.5,'MarkerSize',8);
hold on
plot(SNR_dB,NMSE_dB_1D,'-v', 'Color', colors(3,:),'LineWidth',1.5,'MarkerSize',8);
hold on
plot(SNR_dB,NMSE_dB_4D,'-*', 'Color', colors(2,:),'LineWidth',1.5,'MarkerSize',8);
hold on
plot(SNR_dB,NMSE_dB_omp,'-o', 'Color', colors(5,:),'LineWidth',1.5,'MarkerSize',8);
hold on
plot(SNR_dB,NMSE_dB_4Dsomp,'-^', 'Color', colors(6,:),'LineWidth',1.5,'MarkerSize',8);
grid on

legend('Traditional impulse','1D SBL', '4D SBL','OMP','4D-SOMP');
xlabel('pilot quantity')
ylabel('NMSE(dB)')

%% SNR
% 数据定义
SNR_dB = [0, 5, 10, 15, 20];

% 各算法的NMSE结果 (dB)
NMSE_1D_SBL = [-4.60, -14.52, -22.82, -27.29, -27.66];
NMSE_2D_SBL = [-7.12, -15.57, -19.98, -22.07, -22.45];
NMSE_Hierarchical_1D = [-7.13, -17.18, -26.25, -30.52, -31.36];
NMSE_Hierarchical_2D = [-7.49, -15.94, -20.96, -23.34, -23.86];
NMSE_OMP = [-5.06, -14.34, -19.67, -21.39, -21.46];
NMSE_Traditional = [-4.42, -11.76, -14.08, -14.42, -14.46];
NMSE_CRLB = [-23.79, -29.45, -33.98, -38.73, -43.17];

% 颜色定义
colors = [0,107,182;     % 蓝色 - 1D SBL
          118,41,133;    % 紫色 - 2D SBL  
          234,174,31;    % 黄色 - 分层1D SBL
          215,94,59;     % 橙色 - 分层2D SBL
          184,125,182;   % 淡紫色 - OMP
          71,90,40;      % 绿色 - Traditional
          161,27,30]/256; % 红色 - CRLB

% 创建图形
figure('Position', [100, 100, 800, 600]);
hold on; grid on; box on;

% 绘制各算法性能曲线
plot(SNR_dB, NMSE_1D_SBL, '-s', 'Color', colors(1,:), 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', '1D SBL');
plot(SNR_dB, NMSE_2D_SBL, '-o', 'Color', colors(2,:), 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', '4D SBL');
plot(SNR_dB, NMSE_Hierarchical_1D, '-^', 'Color', colors(3,:), 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Hierarchical 1D SBL');
plot(SNR_dB, NMSE_Hierarchical_2D, '-v', 'Color', colors(4,:), 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Hierarchical 4D SBL');
plot(SNR_dB, NMSE_OMP, '-d', 'Color', colors(5,:), 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'OMP');
plot(SNR_dB, NMSE_Traditional, '-p', 'Color', colors(6,:), 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Traditional Impulse');
plot(SNR_dB, NMSE_CRLB, '--', 'Color', colors(7,:), 'LineWidth', 2, 'DisplayName', 'CRLB');

% 图形设置
xlabel('SNR (dB)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('NMSE (dB)', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 10);

%% speed
% 数据定义
k_max_values = [2, 4, 6, 8, 10];
% 各算法的NMSE结果 (dB)
NMSE_1D_SBL = [-24.17, -22.88, -22.69, -23.83, -22.99];
NMSE_2D_SBL = [-18.96, -19.85, -19.23, -21.31, -20.88];
NMSE_Hierarchical_1D = [-24.91, -24.20, -24.98, -25.40, -24.21];
NMSE_Hierarchical_2D = [-19.43, -19.81, -19.83, -21.25, -21.60];
NMSE_OMP = [-19.82, -19.14, -18.94, -19.75, -19.33];
NMSE_Traditional = [-13.02, -14.02, -14.09, -15.53, -15.11];

% 颜色定义
colors = [0,107,182;     % 蓝色 - 1D SBL
          118,41,133;    % 紫色 - 2D SBL
          234,174,31;    % 黄色 - 分层1D SBL
          215,94,59;     % 橙色 - 分层2D SBL
          184,125,182;   % 淡紫色 - OMP
          71,90,40;      % 绿色 - Traditional
          161,27,30]/256; % 红色 - CRLB

% 创建图形
figure('Position', [100, 100, 800, 600]);
hold on; grid on; box on;

% 绘制各算法性能曲线
plot(k_max_values, NMSE_1D_SBL, '-s', 'Color', colors(1,:), 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', '1D SBL');
plot(k_max_values, NMSE_2D_SBL, '-o', 'Color', colors(2,:), 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', '2D SBL');
plot(k_max_values, NMSE_Hierarchical_1D, '-^', 'Color', colors(3,:), 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', '分层1D SBL');
plot(k_max_values, NMSE_Hierarchical_2D, '-v', 'Color', colors(4,:), 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', '分层2D SBL');
plot(k_max_values, NMSE_OMP, '-d', 'Color', colors(5,:), 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'OMP');
plot(k_max_values, NMSE_Traditional, '-p', 'Color', colors(6,:), 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Traditional Impulse');

% 图形设置
xlabel('k_{max}', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('NMSE (dB)', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 10);
% xlim([-26, -10]);
