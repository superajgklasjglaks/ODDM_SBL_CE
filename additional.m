figure()
% hold on;
% plot(SNR_dB,NMSE_dB_trd,'-s', 'Color', colors(1,:),'LineWidth',1.5,'MarkerSize',8);
% plot(SNR_dB,NMSE_dB_1D,'-v', 'Color', colors(3,:),'LineWidth',1.5,'MarkerSize',8);
% plot(SNR_dB,NMSE_dB_2D,'-*', 'Color', colors(2,:),'LineWidth',1.5,'MarkerSize',8);
% plot(SNR_dB,NMSE_dB_omp,'-o', 'Color', colors(6,:),'LineWidth',1.5,'MarkerSize',8);
% plot(SNR_dB,NMSE_dB_CRLB,'-^', 'Color', colors(5,:), 'LineWidth',1.5,'MarkerSize',8);
% % hold on
% % plot(SNR_dB,NMSE_dB_1D_o,'-o', 'Color', colors(6,:),'LineWidth',1.5,'MarkerSize',8);
% % hold on
% % plot(SNR_dB,NMSE_dB_2D_o,'-^', 'Color', colors(7,:),'LineWidth',1.5,'MarkerSize',8);
% grid on
% 
% legend('Traditional impulse','1D off-grid','2D off-grid','OMP','CRLB');
% xlabel('SNR(dB)')
% ylabel('NMSE(dB)')

hold on;
plot(SNR_dB, NMSE_dB_trd, '-s', 'Color', colors(1,:), 'LineWidth', 1.5, 'MarkerSize', 8, 'DisplayName', 'Traditional');
plot(SNR_dB, NMSE_dB_1D, '-v', 'Color', colors(2,:), 'LineWidth', 1.5, 'MarkerSize', 8, 'DisplayName', '1D');
plot(SNR_dB, NMSE_dB_2D, '-*', 'Color', colors(3,:), 'LineWidth', 1.5, 'MarkerSize', 8, 'DisplayName', '2D');
plot(SNR_dB, NMSE_dB_omp, '-o', 'Color', colors(4,:), 'LineWidth', 1.5, 'MarkerSize', 8, 'DisplayName', 'OMP');
plot(SNR_dB, NMSE_dB_CRLB, '-^', 'Color', colors(5,:), 'LineWidth', 1.5, 'MarkerSize', 8, 'DisplayName', 'CRLB');
grid on;

legend('Location', 'best'); % 自动匹配 DisplayName
xlabel('SNR (dB)');
ylabel('NMSE (dB)');