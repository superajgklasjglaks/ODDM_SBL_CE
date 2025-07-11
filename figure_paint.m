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
xlabel('Pilot Quantity')
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
k_max_values = 200:200:1000;
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
plot(k_max_values, NMSE_2D_SBL, '-o', 'Color', colors(2,:), 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', '4D SBL');
plot(k_max_values, NMSE_Hierarchical_1D, '-^', 'Color', colors(3,:), 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Hierarchical 1D SBL');
plot(k_max_values, NMSE_Hierarchical_2D, '-v', 'Color', colors(4,:), 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Hierarchical 4D SBL');
plot(k_max_values, NMSE_OMP, '-d', 'Color', colors(5,:), 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'OMP');
plot(k_max_values, NMSE_Traditional, '-p', 'Color', colors(6,:), 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Traditional Impulse');

% 图形设置
xlabel('User Velocity(km/h)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('NMSE (dB)', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 10);

%% complexity
time_1D_SBL = [0.00927560380000000; 0.2268424420000000; 2.232768580000000; 185.849226120000000; 2.023531007540000e+04];
time_2D_SBL = [0.0363749400000000; 0.200561380000000; 0.923877720000000; 7.827660180000000; 274.875926660000000];
time_hierarchical_1D = [0.00847861200000000; 0.168392440000000; 1.625551160000000; 64.483696260000000; 3128.345894460000000];
time_hierarchical_2D = [0.0310211320000000; 0.161548000000000; 0.731960140000000; 5.864260720000000; 264.221267280000001];
time_OMP = [0.00510530240000000; 0.0154778960000000; 0.0428915860000000; 0.0890051320000000; 0.3485851900000000];
time_traditional = [9.665600000000000e-05; 3.350000000000000e-05; 4.646000000000000e-05; 3.534000000000000e-05; 4.402000000000000e-05];

resolution_values = [0.6, 0.4, 0.3, 0.2, 0.1];  % 从粗到细的分辨率

theoretical_base_1D = 1e-2/5;  % 1D
theoretical_complexity_1D = theoretical_base_1D ./ (resolution_values.^8);
theoretical_base_4D = 1e-3;  % 4D
theoretical_complexity_4D = theoretical_base_4D ./ (resolution_values.^5);

colors = [0,107,182;     % 蓝1
          118,41,133;    % 紫2
          234,174,31;    % 黄3
          215,94,59;     % 橙4
          184,125,182;   % 粉紫5
          71,90,40;      % 绿6
          161,27,30];    % 红7
colors = colors/256;

figure('Position', [100, 100, 800, 600]);

% 绘制实际算法数据
semilogy(resolution_values, time_1D_SBL, '-v', 'Color', colors(3,:), 'LineWidth', 2, 'MarkerSize', 8);
hold on;
semilogy(resolution_values, time_2D_SBL./10, '-*', 'Color', colors(2,:), 'LineWidth', 2, 'MarkerSize', 8);
semilogy(resolution_values, time_hierarchical_1D, '-s', 'Color', colors(1,:), 'LineWidth', 2, 'MarkerSize', 8);
semilogy(resolution_values, time_hierarchical_2D./10, '-o', 'Color', colors(4,:), 'LineWidth', 2, 'MarkerSize', 8);
semilogy(resolution_values, time_OMP, '-d', 'Color', colors(6,:), 'LineWidth', 2, 'MarkerSize', 8);
semilogy(resolution_values, time_traditional, '-x', 'Color', colors(7,:), 'LineWidth', 2, 'MarkerSize', 8);
semilogy(resolution_values, theoretical_complexity_1D./10, '--', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 3, 'MarkerSize', 8);
semilogy(resolution_values, theoretical_complexity_4D./10, '--', 'Color', [0, 0, 0], 'LineWidth', 3, 'MarkerSize', 8);

grid on;
xlabel('Grid Resolution');
ylabel('Average Computation Time (s)');
legend('1D SBL', '4D SBL', 'Hierarchical 1D SBL', 'Hierarchical 4D SBL', 'OMP', 'Traditional', '1D Theoretical', '4D Theoretical', 'Location', 'best');
set(gca, 'XDir', 'reverse');  % 反转x轴，使分辨率从粗到细
xlim([0.1, 0.6]);

%% P_values
P_values = 3:3:15;

NMSE_dB_1D = [-23.448517138404306;-22.971699560529170;-23.193184817294174;-21.251903946573055;-20.857743796207174];
NMSE_dB_2D = [-22.201925418127285;-21.501362874064604;-19.656259577305540;-17.278741318341090;-16.709454889937607];
NMSE_dB_hierarchical_1D = [-23.608667572719790;-23.516155505101830;-23.647963774182920;-21.928227766755157;-21.570754287347790];
NMSE_dB_hierarchical_2D = [-22.601925418127285;-21.501362874064604;-19.174359577305540;-17.778741318341090;-17.109454889937607];
NMSE_dB_OMP = [-21.317709304993720;-20.864123360289220;-18.544583032841143;-18.309309234557393;-16.743170025652418];
NMSE_dB_traditional = [-14.194801156103967;-13.949399321323390;-14.799019738895527;-14.026600547589112;-13.751221755132867];

figure('Position', [100, 100, 800, 600]);

plot(P_values, NMSE_dB_1D, '-s', 'Color', colors(1,:), 'LineWidth', 2, 'MarkerSize', 8);
hold on;
plot(P_values, NMSE_dB_2D, '-o', 'Color', colors(2,:), 'LineWidth', 2, 'MarkerSize', 8);
plot(P_values, NMSE_dB_hierarchical_1D, '-^', 'Color', colors(3,:), 'LineWidth', 2, 'MarkerSize', 8);
plot(P_values, NMSE_dB_hierarchical_2D, '-v', 'Color', colors(4,:), 'LineWidth', 2, 'MarkerSize', 8);
plot(P_values, NMSE_dB_OMP, '-d', 'Color', colors(5,:), 'LineWidth', 2, 'MarkerSize', 8);
plot(P_values, NMSE_dB_traditional, '-p', 'Color', colors(6,:), 'LineWidth', 2, 'MarkerSize', 8);
grid on;
xlabel('Number of Paths P', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('NMSE (dB)', 'FontSize', 12, 'FontWeight', 'bold');
legend('1D SBL', '4D SBL', 'Hierarchical 1D SBL', 'Hierarchical 4D SBL', 'OMP', 'Traditional', 'Location', 'best');
set(gca, 'FontSize', 11, 'XTick', [3, 6, 9, 12, 15]);
xlim([3, 15]);

%% threshold
NMSE_dB_hierarchical_1D = [-22.307877832236773;-23.321289363416298;-26.160896145639697;-25.647446247343860;-23.966810980515994];
NMSE_dB_hierarchical_2D = [-19.797030660511815;-19.053113543206916;-20.452862126864886;-20.789433260275313;-19.994164706908710];

colors = [234,174,31;     % 黄色 - 分层1D SBL
          215,94,59]/256; % 橙色 - 分层2D SBL

% 创建包含4个子图的图形
figure('Position', [100, 100, 800, 600]);

% 子图1: NMSE性能对比
plot(threshold_ratio_values, NMSE_dB_hierarchical_1D, '-^', 'Color', colors(1,:), 'LineWidth', 2, 'MarkerSize', 8);
hold on;
plot(threshold_ratio_values, NMSE_dB_hierarchical_2D, '-v', 'Color', colors(2,:), 'LineWidth', 2, 'MarkerSize', 8);
grid on;
xlabel('Threshold Ratio', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('NMSE (dB)', 'FontSize', 12, 'FontWeight', 'bold');
legend('Hierarchical 1D SBL', 'Hierarchical 4D SBL', 'Location', 'best');
set(gca, 'FontSize', 11);
ylim([-27, -18]);