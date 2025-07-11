% ����Massive MIMO�Ƕȹ��ƴ���
close all
clear all
clc

fprintf('��ʼ����Massive MIMO�Ƕȹ��ƴ���...\n');

try
    % ���ü򻯵Ĳ��Բ���
    Nt = 8;  % �����������Լӿ����
    Nr = 8;
    P = 2;   % ����·����
    N_fram = 2;  % ���ٷ���֡��
    
    % �ǶȲ���
    aod_max = pi/6;  % 30��
    aoa_max = pi/6;
    
    % �ֱ��ʲ���
    r_aod = 0.05;
    r_aoa = 0.05;
    
    % SNR����
    SNR_dB = [0, 10];
    SNR = 10.^(SNR_dB/10);
    sigma_2 = 0.5 ./ SNR;
    
    fprintf('���Բ����������\n');
    
    % ����������Ӧ����
    fprintf('����������Ӧ����...\n');
    aod_test = 0.1;
    aoa_test = 0.2;
    
    w_tx = MIMO_Array_Response_Tx(Nt, 1, 1, aod_test);
    w_rx = MIMO_Array_Response_Rx(Nr, 1, 1, aoa_test);
    
    fprintf('����������Ӧ: %.4f + %.4fi\n', real(w_tx), imag(w_tx));
    fprintf('����������Ӧ: %.4f + %.4fi\n', real(w_rx), imag(w_rx));
    
    % ���ԽǶ���������
    fprintf('���ԽǶ���������...\n');
    N_aod = ceil(2 * aod_max / r_aod);
    M_aoa = ceil(2 * aoa_max / r_aoa);
    
    [aod_bar, aoa_bar] = MIMO_First_Order_Linear_Approximation(N_aod, M_aoa, aod_max, r_aod, r_aoa);
    fprintf('���ɽǶ�����: %d x %d = %d ����\n', N_aod, M_aoa, length(aod_bar));
    
    % ���Լ򵥵�1D SBL
    fprintf('����1D SBL�㷨...\n');
    
    % ���ɲ�������
    aod_true = [0.1, -0.15];
    aoa_true = [0.2, -0.1];
    h_true = [1, 0.5];
    
    % ������������
    M_T = 4;
    N_T = 4;
    y_T = zeros(M_T * N_T, 1);
    
    for p = 1:length(aod_true)
        for nnn = 0:(N_T-1)
            for mmm = 1:M_T
                matrix_idx = nnn * M_T + mmm;
                y_T(matrix_idx) = y_T(matrix_idx) + h_true(p) * ...
                    MIMO_Array_Response_Tx(Nt, nnn, 3, aod_true(p)) * ...
                    MIMO_Array_Response_Rx(Nr, mmm-1, 1, aoa_true(p));
            end
        end
    end
    
    % �������
    noise = (sigma_2(1) * randn(size(y_T)) + 1j * sigma_2(1) * randn(size(y_T)));
    y_T = y_T + noise;
    
    % ����1D SBL
    [h_hat, aod_hat, aoa_hat, virtual_size, phi_trunc, delta_a] = ...
        MIMO_CE_1D_SBL(1, 4, 2, Nr, Nt, N_T, M_T, y_T, r_aod, r_aoa, aod_max, aoa_max, 0);
    
    fprintf('1D SBL��ɣ������� %d ��ϵ��\n', length(h_hat));
    fprintf('���ϵ������: %.4f\n', max(abs(h_hat)));
    
    % ����OMP
    fprintf('����OMP�㷨...\n');
    [h_hat_omp, omp_index, aod_hat_omp, aoa_hat_omp] = ...
        MIMO_OMP(1, 4, 2, Nr, Nt, N_T, M_T, y_T, r_aod, r_aoa, aod_max, aoa_max);
    
    fprintf('OMP��ɣ�ѡ���� %d ��ԭ��\n', length(omp_index));
    
    % ���Դ�ͳ��������
    fprintf('���Դ�ͳ��������...\n');
    y_trunc = reshape(y_T, M_T, N_T).';
    [h_hat_trad, aod_hat_trad, aoa_hat_trad] = ...
        MIMO_traditional_beamforming(1, y_trunc, 2, 2, sigma_2(1));
    
    fprintf('��ͳ����������ɣ������� %d ��ϵ��\n', length(h_hat_trad));
    
    fprintf('\n=== ���в���ͨ��! ===\n');
    fprintf('��������������У�����ִ�������ķ��档\n');
    
catch ME
    fprintf('\n=== ����ʧ�� ===\n');
    fprintf('������Ϣ: %s\n', ME.message);
    fprintf('����λ��: %s (�� %d ��)\n', ME.stack(1).name, ME.stack(1).line);
    
    % ��ʾ��ϸ�Ĵ����ջ
    for i = 1:length(ME.stack)
        fprintf('  %d. %s (�� %d ��)\n', i, ME.stack(i).name, ME.stack(i).line);
    end
end