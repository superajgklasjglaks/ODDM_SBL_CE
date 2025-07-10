% clc;
% clear;
% %% 初始化及设定参数
% array_num = 10;%%阵元数
% snapshot_num = 100;%%快拍数
% source_aoa = [-30,-5,45];%%信源到达角
% c = 340;%%波速
% f = 1000;%%频率
% lambda = c/f;%%波长
% d = 0.5*lambda;
% source_num = length(source_aoa);%%信源数
% sig_nr = [20,20,20];%%信噪比、扰噪比
% reso_num = 91;%%网格数
% %% 导向矢量
% X = zeros(source_num,snapshot_num);
% A = exp(-1i*(0:array_num-1)'*2*pi*(d/lambda)*sind(source_aoa));%%阵列响应矩阵
% for ik = 1:length(sig_nr)
%      X(ik,:) = sqrt(10^(sig_nr(ik)/10))*(randn(1,snapshot_num)+randn(1,snapshot_num)*1i)/sqrt(2);
% end
% n = (randn(array_num,snapshot_num)+randn(array_num,snapshot_num)*1i)/sqrt(2);
% Y = A*X+n;
% % [~,~,D_svd] = svd(Y,'econ');
% % Y = Y*D_svd(:,1:source_num);%%信号降维
% %% OGSBI算法输入量整理
% params.Y = Y;
% params.reso_num = reso_num;
% params.maxiter = 2000;%%最大迭代次数
% params.tolerance = 1e-4;%%误差容忍度
% params.sigma2 = mean(var(Y))/100;%%噪声方差估计值
% res = OGSBI(params);
% xp_rec = res.reso_grid;
% x_rec = res.mu;
% % x_rec = res.mu * D_svd(:,1:size(res.mu,2))';
% xpower_rec = mean(abs(x_rec).^2,2) + real(diag(res.Sigma)) * source_num / snapshot_num;
% xpower_rec = abs(xpower_rec)/max(xpower_rec);
% [xp_rec,x_index] = sort(xp_rec,'ascend');
% xpower_rec = xpower_rec(x_index);
% figure();
% plot(xp_rec,10*log(xpower_rec));xlabel("角度/°");ylabel("归一化功率/dB");
% hold on;
% semilogy(source_aoa,max(xpower_rec),'bo');
% hold off;
% 
% function res = OGSBI(params)
% %% 函数参数初始化
% Y = params.Y;
% reso_num = params.reso_num;
% reso_grid = linspace(-90,90,reso_num)';
% reso = 180/(reso_num-1);
% [array_num, snapshot_num] = size(Y);
% r = reso*pi/180;
% maxiter = params.maxiter;
% tol = params.tolerance;
% index_b = randperm(length(reso_grid),array_num)';%%该变量主要记录alpha中最大的几个元素的位置，以后续对这几个位置进行一阶泰勒展开
% converged = false;%%判断收敛的Boolen
% iter = 0;
% A = exp(-1i*(0:array_num-1)'*pi*sind(reso_grid'));
% B = -1i*pi*(0:array_num-1)'*cosd(reso_grid').*A;
% alpha = mean(abs(A'*Y), 2);
% beta = zeros(reso_num,1);
% c_sigma0_init = 1e-4;
% d_sigma0_init = 1e-4;
% c_gamma_init = 1;
% d_gamma_init = 1e-4;
% alpha_0 = 0.01;
% alpha_0_seq = zeros(maxiter,1);%%噪声精度变化迭代
% Phi = A;
% while ~converged
%     iter = iter+1;
%     Phi(:,index_b) = exp(-1i*(0:array_num-1)'*pi*sind(reso_grid(index_b)'));%%根据上一轮迭代得到的网格点位置，对这几个点进行进一步迭代。
%     B(:,index_b) = -1i*pi*(0:array_num-1)'*cosd(reso_grid(index_b)').*A(:,index_b);
%     alpha_last = alpha;%%上一次迭代出的alpha的结果
% %% 更新X的后验概率密度函数的均值与方差
%     C = 1/alpha_0*eye(array_num)+Phi*diag(alpha)*Phi';
%     Cinv = inv(C);%%(16)的woodbury反演形式中的逆矩阵的括号内部分
%     Sigma = diag(alpha)-diag(alpha)*Phi'*Cinv*Phi*diag(alpha);%%(16)的woodbury矩阵反演公式形式
%     mu = alpha_0*Sigma*Phi'*Y;%%(15)
% %% 更新alpha
%     musq = mean(abs(mu).^2,2);
%     alpha = musq+real(diag(Sigma));
%     for ik = 1:reso_num
%         alpha(ik) = (-snapshot_num+sqrt(snapshot_num^2+4*d_gamma_init*(mu(ik,:)*mu(ik,:)'+snapshot_num*Sigma(ik,ik))))/(2*d_gamma_init);
%     end
% %% 更新alpha_0
%     alpha_0 = (snapshot_num*array_num+c_sigma0_init-1)/(norm(Y-Phi*mu,'fro')^2+snapshot_num*trace(Phi*Sigma*Phi')+d_sigma0_init);%%(18),范数部分详见18至19之间部分
%     alpha_0_seq(iter) = alpha_0;
% %% 判断是否停止迭代
%     if norm(alpha-alpha_last)/norm(alpha_last) < tol || iter >= maxiter
%         converged = true;
%     end%%收敛或迭代次数达到上限时进入该循环
% %% 更新beta
%     [beta,index_b] = off_grid_operation(Y,alpha,array_num,mu,Sigma,Phi,B,beta,r);
%     reso_grid(index_b) = reso_grid(index_b)+beta(index_b)*180/pi;
% end
% res.mu = mu;
% res.Sigma = Sigma;
% res.beta = beta;
% res.alpha = alpha;
% res.iter = iter;
% res.sigma2 = 1/alpha_0;
% res.sigma2seq = 1./alpha_0_seq(1:iter);
% res.reso_grid = reso_grid';
% end
% 
% function [beta,index_b] = off_grid_operation(Y,gamma,iter_size,mu,Sigma,Phi,B,beta,r)
%     reso_num = size(B,2);
%     snapshot_num = size(Y,2);
%     [~, index_b] = sort(gamma, 'descend');
%     index_b = index_b(1:iter_size);%%选定位置进行一阶泰勒展开
%     temp = beta;
%     beta = zeros(reso_num,1);
%     beta(index_b) = temp(index_b);
%     BHB = B'*B;
%     P = real(conj(BHB(index_b,index_b)).*(mu(index_b,:)*mu(index_b,:)'+snapshot_num*Sigma(index_b,index_b)));%%(20)
%     v = zeros(length(index_b), 1);%%(21)
%     for t = 1:snapshot_num
%         v = v+real(conj(mu(index_b,t)).*(B(:,index_b)'*(Y(:,t)-Phi*mu(:,t))));
%     end
%     v = v-snapshot_num*real(diag(B(:,index_b)'*Phi*Sigma(:,index_b)));
%     eigP = svd(P);
%     eigP = sort(eigP,'descend');
%     if eigP(end)/eigP(1) > 1e-5 || any(diag(P) == 0)
%         for n = 1:iter_size
%             temp_beta = beta(index_b);
%             temp_beta(n) = 0;
%             beta(index_b(n)) = (v(n)-P(n,:)*temp_beta)/P(n,n);%%(26.1)
%             if abs(beta(index_b(n))) > r/2%%(26.2)
%                 beta(index_b(n)) = r/2*sign(beta(index_b(n)));
%             end
%             if P(n,n) == 0
%                 beta(index_b(n)) = 0;
%             end
%         end
%     else
%         beta = zeros(reso_num,1);
%         beta(index_b) = P\v;
%     end  
% end

clc;
clear;
%% 初始化及设定参数
array_num = 10;              % 阵元数
snapshot_num = 100;          % 快拍数
source_aoa = [-30, -5, 10, 45];  % 信源真实角度
c = 340;                     % 波速
f = 1000;                    % 频率
lambda = c / f;              % 波长
d = 0.5 * lambda;            % 阵元间距
source_num = length(source_aoa); % 信源数
reso_num = 31;               % 网格数
SNR_range = -10:5:20;        % 信噪比范围（dB）
MC_trials = 50;              % 蒙特卡洛仿真次数

%% 存储NMSE结果
NMSE = zeros(length(SNR_range), 1);

%% 主循环：不同SNR
for snr_idx = 1:length(SNR_range)
    SNR = SNR_range(snr_idx);
    fprintf('Processing SNR = %d dB...\n', SNR);

    % 临时存储单次SNR下的误差
    trial_errors = zeros(MC_trials, 1);

    %% 蒙特卡洛仿真
    for mc = 1:MC_trials
        %% 生成接收信号
        A = exp(-1i * (0:array_num-1)' * 2 * pi * (d / lambda) * sind(source_aoa)); % 阵列响应矩阵
        X = (randn(source_num, snapshot_num) + 1i * randn(source_num, snapshot_num)) / sqrt(2); % 信源信号
        signal_power = mean(abs(A * X).^2, 'all'); % 信号功率
        noise_power = signal_power / (10^(SNR / 10)); % 噪声功率
        n = sqrt(noise_power) * (randn(array_num, snapshot_num) + 1i * randn(array_num, snapshot_num)) / sqrt(2);
        Y = A * X + n; % 接收信号

        %% OGSBI算法
        params.Y = Y;
        params.reso_num = reso_num;
        params.maxiter = 2000;
        params.tolerance = 1e-4;
        params.sigma2 = noise_power; % 使用真实噪声功率初始化
        res = OGSBI(params);

        %% 提取估计角度
        xpower_rec = mean(abs(res.mu).^2, 2) + real(diag(res.Sigma)) * source_num / snapshot_num;
        xpower_rec = abs(xpower_rec) / max(xpower_rec);
        [~, est_peaks] = findpeaks(xpower_rec, 'SortStr', 'descend', 'NPeaks', source_num);
        est_aoa = res.reso_grid(est_peaks);

        %% 计算角度误差（匹配最近的真值）
        error = 0;
        for k = 1:source_num
            [~, idx] = min(abs(est_aoa - source_aoa(k)));
            error = error + (est_aoa(idx) - source_aoa(k))^2;
            est_aoa(idx) = Inf; % 避免重复匹配
        end
        trial_errors(mc) = error / source_num;
    end

    %% 计算NMSE
    NMSE(snr_idx) = mean(trial_errors) / (mean(source_aoa.^2)); % 归一化
end

%% 绘制NMSE vs SNR
figure;
plot(SNR_range, 10*log10(NMSE), 'bo-', 'LineWidth', 2);
grid on;
xlabel('SNR (dB)');
ylabel('NMSE');
title('DOA Estimation NMSE vs SNR');

function res = OGSBI(params)
%% 函数参数初始化
Y = params.Y;
reso_num = params.reso_num;
reso_grid = linspace(-90,90,reso_num)';
reso = 180/(reso_num-1);
[array_num, snapshot_num] = size(Y);
r = reso*pi/180;
maxiter = params.maxiter;
tol = params.tolerance;
index_b = randperm(length(reso_grid),array_num)';%%该变量主要记录alpha中最大的几个元素的位置，以后续对这几个位置进行一阶泰勒展开
converged = false;%%判断收敛的Boolen
iter = 0;
A = exp(-1i*(0:array_num-1)'*pi*sind(reso_grid'));
B = -1i*pi*(0:array_num-1)'*cosd(reso_grid').*A;
alpha = mean(abs(A'*Y), 2);
beta = zeros(reso_num,1);
c_sigma0_init = 1e-4;
d_sigma0_init = 1e-4;
c_gamma_init = 1;
d_gamma_init = 1e-4;
alpha_0 = 0.01;
alpha_0_seq = zeros(maxiter,1);%%噪声精度变化迭代
Phi = A;
while ~converged
    iter = iter+1;
    Phi(:,index_b) = exp(-1i*(0:array_num-1)'*pi*sind(reso_grid(index_b)'));%%根据上一轮迭代得到的网格点位置，对这几个点进行进一步迭代。
    B(:,index_b) = -1i*pi*(0:array_num-1)'*cosd(reso_grid(index_b)').*A(:,index_b);
    alpha_last = alpha;%%上一次迭代出的alpha的结果
%% 更新X的后验概率密度函数的均值与方差
    C = 1/alpha_0*eye(array_num)+Phi*diag(alpha)*Phi';
    Cinv = inv(C);%%(16)的woodbury反演形式中的逆矩阵的括号内部分
    Sigma = diag(alpha)-diag(alpha)*Phi'*Cinv*Phi*diag(alpha);%%(16)的woodbury矩阵反演公式形式
    mu = alpha_0*Sigma*Phi'*Y;%%(15)
%% 更新alpha
    musq = mean(abs(mu).^2,2);
    alpha = musq+real(diag(Sigma));
    for ik = 1:reso_num
        alpha(ik) = (-snapshot_num+sqrt(snapshot_num^2+4*d_gamma_init*(mu(ik,:)*mu(ik,:)'+snapshot_num*Sigma(ik,ik))))/(2*d_gamma_init);
    end
%% 更新alpha_0
    alpha_0 = (snapshot_num*array_num+c_sigma0_init-1)/(norm(Y-Phi*mu,'fro')^2+snapshot_num*trace(Phi*Sigma*Phi')+d_sigma0_init);%%(18),范数部分详见18至19之间部分
    alpha_0_seq(iter) = alpha_0;
%% 判断是否停止迭代
    if norm(alpha-alpha_last)/norm(alpha_last) < tol || iter >= maxiter
        converged = true;
    end%%收敛或迭代次数达到上限时进入该循环
%% 更新beta
    [beta,index_b] = off_grid_operation(Y,alpha,array_num,mu,Sigma,Phi,B,beta,r);
    reso_grid(index_b) = reso_grid(index_b)+beta(index_b)*180/pi;
end
res.mu = mu;
res.Sigma = Sigma;
res.beta = beta;
res.alpha = alpha;
res.iter = iter;
res.sigma2 = 1/alpha_0;
res.sigma2seq = 1./alpha_0_seq(1:iter);
res.reso_grid = reso_grid';
end

function [beta,index_b] = off_grid_operation(Y,gamma,iter_size,mu,Sigma,Phi,B,beta,r)
    reso_num = size(B,2);
    snapshot_num = size(Y,2);
    [~, index_b] = sort(gamma, 'descend');
    index_b = index_b(1:iter_size);%%选定位置进行一阶泰勒展开
    temp = beta;
    beta = zeros(reso_num,1);
    beta(index_b) = temp(index_b);
    BHB = B'*B;
    P = real(conj(BHB(index_b,index_b)).*(mu(index_b,:)*mu(index_b,:)'+snapshot_num*Sigma(index_b,index_b)));%%(20)
    v = zeros(length(index_b), 1);%%(21)
    for t = 1:snapshot_num
        v = v+real(conj(mu(index_b,t)).*(B(:,index_b)'*(Y(:,t)-Phi*mu(:,t))));
    end
    v = v-snapshot_num*real(diag(B(:,index_b)'*Phi*Sigma(:,index_b)));
    eigP = svd(P);
    eigP = sort(eigP,'descend');
    if eigP(end)/eigP(1) > 1e-5 || any(diag(P) == 0)
        for n = 1:iter_size
            temp_beta = beta(index_b);
            temp_beta(n) = 0;
            beta(index_b(n)) = (v(n)-P(n,:)*temp_beta)/P(n,n);%%(26.1)
            if abs(beta(index_b(n))) > r/2%%(26.2)
                beta(index_b(n)) = r/2*sign(beta(index_b(n)));
            end
            if P(n,n) == 0
                beta(index_b(n)) = 0;
            end
        end
    else
        beta = zeros(reso_num,1);
        beta(index_b) = P\v;
    end  
end

% clc;
% clear;
% %% 初始化及设定参数
% array_nums = 8:8:32;         % 阵元数范围
% snapshot_num = 100;          % 快拍数
% source_aoa = [-30, -5, 10, 45];  % 信源真实角度
% c = 340;                     % 波速
% f = 1000;                    % 频率
% lambda = c / f;              % 波长
% d = 0.5 * lambda;            % 阵元间距
% source_num = length(source_aoa); % 信源数
% SNR = 10;                    % 固定信噪比为10dB
% MC_trials = 50;              % 蒙特卡洛仿真次数
% 
% %% 存储NMSE结果
% NMSE = zeros(length(array_nums), 1);
% 
% %% 主循环：不同天线数
% for array_idx = 1:length(array_nums)
%     array_num = array_nums(array_idx);
%     fprintf('Processing array number = %d...\n', array_num);
% 
%     % 动态设置网格数（至少为天线数的3倍）
%     reso_num = max(31, 3*array_num);  % 保证最小网格数为31，且至少为天线数的3倍
% 
%     % 临时存储单次天线数下的误差
%     trial_errors = zeros(MC_trials, 1);
% 
%     %% 蒙特卡洛仿真
%     for mc = 1:MC_trials
%         %% 生成接收信号
%         A = exp(-1i * (0:array_num-1)' * 2 * pi * (d / lambda) * sind(source_aoa)); % 阵列响应矩阵
%         X = (randn(source_num, snapshot_num) + 1i * randn(source_num, snapshot_num)) / sqrt(2); % 信源信号
%         signal_power = mean(abs(A * X).^2, 'all'); % 信号功率
%         noise_power = signal_power / (10^(SNR / 10)); % 噪声功率
%         n = sqrt(noise_power) * (randn(array_num, snapshot_num) + 1i * randn(array_num, snapshot_num)) / sqrt(2);
%         Y = A * X + n; % 接收信号
% 
%         %% OGSBI算法
%         params.Y = Y;
%         params.reso_num = reso_num;
%         params.maxiter = 2000;
%         params.tolerance = 1e-4;
%         params.sigma2 = noise_power; % 使用真实噪声功率初始化
%         res = OGSBI(params);
% 
%         %% 提取估计角度
%         xpower_rec = mean(abs(res.mu).^2, 2) + real(diag(res.Sigma)) * source_num / snapshot_num;
%         xpower_rec = abs(xpower_rec) / max(xpower_rec);
%         [~, est_peaks] = findpeaks(xpower_rec, 'SortStr', 'descend', 'NPeaks', source_num);
%         est_aoa = res.reso_grid(est_peaks);
% 
%         %% 计算角度误差（匹配最近的真值）
%         error = 0;
%         for k = 1:source_num
%             [~, idx] = min(abs(est_aoa - source_aoa(k)));
%             error = error + (est_aoa(idx) - source_aoa(k))^2;
%             est_aoa(idx) = Inf; % 避免重复匹配
%         end
%         trial_errors(mc) = error / source_num;
%     end
% 
%     %% 计算NMSE
%     NMSE(array_idx) = mean(trial_errors) / (mean(source_aoa.^2)); % 归一化
% end
% 
% %% 绘制NMSE vs 天线数
% figure;
% plot(array_nums, 10*log10(NMSE), 'bo-', 'LineWidth', 2);
% grid on;
% xlabel('Number of Antennas');
% ylabel('NMSE (dB)');
% title('DOA Estimation NMSE vs Number of Antennas at SNR=10dB');
% 
% % OGSBI函数保持不变...
% function res = OGSBI(params)
% %% 函数参数初始化
% Y = params.Y;
% reso_num = params.reso_num;
% reso_grid = linspace(-90,90,reso_num)';
% reso = 180/(reso_num-1);
% [array_num, snapshot_num] = size(Y);
% r = reso*pi/180;
% maxiter = params.maxiter;
% tol = params.tolerance;
% index_b = randperm(length(reso_grid),array_num)';%%该变量主要记录alpha中最大的几个元素的位置，以后续对这几个位置进行一阶泰勒展开
% converged = false;%%判断收敛的Boolen
% iter = 0;
% A = exp(-1i*(0:array_num-1)'*pi*sind(reso_grid'));
% B = -1i*pi*(0:array_num-1)'*cosd(reso_grid').*A;
% alpha = mean(abs(A'*Y), 2);
% beta = zeros(reso_num,1);
% c_sigma0_init = 1e-4;
% d_sigma0_init = 1e-4;
% c_gamma_init = 1;
% d_gamma_init = 1e-4;
% alpha_0 = 0.01;
% alpha_0_seq = zeros(maxiter,1);%%噪声精度变化迭代
% Phi = A;
% while ~converged
%     iter = iter+1;
%     Phi(:,index_b) = exp(-1i*(0:array_num-1)'*pi*sind(reso_grid(index_b)'));%%根据上一轮迭代得到的网格点位置，对这几个点进行进一步迭代。
%     B(:,index_b) = -1i*pi*(0:array_num-1)'*cosd(reso_grid(index_b)').*A(:,index_b);
%     alpha_last = alpha;%%上一次迭代出的alpha的结果
% %% 更新X的后验概率密度函数的均值与方差
%     C = 1/alpha_0*eye(array_num)+Phi*diag(alpha)*Phi';
%     Cinv = inv(C);%%(16)的woodbury反演形式中的逆矩阵的括号内部分
%     Sigma = diag(alpha)-diag(alpha)*Phi'*Cinv*Phi*diag(alpha);%%(16)的woodbury矩阵反演公式形式
%     mu = alpha_0*Sigma*Phi'*Y;%%(15)
% %% 更新alpha
%     musq = mean(abs(mu).^2,2);
%     alpha = musq+real(diag(Sigma));
%     for ik = 1:reso_num
%         alpha(ik) = (-snapshot_num+sqrt(snapshot_num^2+4*d_gamma_init*(mu(ik,:)*mu(ik,:)'+snapshot_num*Sigma(ik,ik))))/(2*d_gamma_init);
%     end
% %% 更新alpha_0
%     alpha_0 = (snapshot_num*array_num+c_sigma0_init-1)/(norm(Y-Phi*mu,'fro')^2+snapshot_num*trace(Phi*Sigma*Phi')+d_sigma0_init);%%(18),范数部分详见18至19之间部分
%     alpha_0_seq(iter) = alpha_0;
% %% 判断是否停止迭代
%     if norm(alpha-alpha_last)/norm(alpha_last) < tol || iter >= maxiter
%         converged = true;
%     end%%收敛或迭代次数达到上限时进入该循环
% %% 更新beta
%     [beta,index_b] = off_grid_operation(Y,alpha,array_num,mu,Sigma,Phi,B,beta,r);
%     reso_grid(index_b) = reso_grid(index_b)+beta(index_b)*180/pi;
% end
% res.mu = mu;
% res.Sigma = Sigma;
% res.beta = beta;
% res.alpha = alpha;
% res.iter = iter;
% res.sigma2 = 1/alpha_0;
% res.sigma2seq = 1./alpha_0_seq(1:iter);
% res.reso_grid = reso_grid';
% end
% 
% function [beta,index_b] = off_grid_operation(Y,gamma,iter_size,mu,Sigma,Phi,B,beta,r)
%     reso_num = size(B,2);
%     snapshot_num = size(Y,2);
%     [~, index_b] = sort(gamma, 'descend');
%     index_b = index_b(1:iter_size);%%选定位置进行一阶泰勒展开
%     temp = beta;
%     beta = zeros(reso_num,1);
%     beta(index_b) = temp(index_b);
%     BHB = B'*B;
%     P = real(conj(BHB(index_b,index_b)).*(mu(index_b,:)*mu(index_b,:)'+snapshot_num*Sigma(index_b,index_b)));%%(20)
%     v = zeros(length(index_b), 1);%%(21)
%     for t = 1:snapshot_num
%         v = v+real(conj(mu(index_b,t)).*(B(:,index_b)'*(Y(:,t)-Phi*mu(:,t))));
%     end
%     v = v-snapshot_num*real(diag(B(:,index_b)'*Phi*Sigma(:,index_b)));
%     eigP = svd(P);
%     eigP = sort(eigP,'descend');
%     if eigP(end)/eigP(1) > 1e-5 || any(diag(P) == 0)
%         for n = 1:iter_size
%             temp_beta = beta(index_b);
%             temp_beta(n) = 0;
%             beta(index_b(n)) = (v(n)-P(n,:)*temp_beta)/P(n,n);%%(26.1)
%             if abs(beta(index_b(n))) > r/2%%(26.2)
%                 beta(index_b(n)) = r/2*sign(beta(index_b(n)));
%             end
%             if P(n,n) == 0
%                 beta(index_b(n)) = 0;
%             end
%         end
%     else
%         beta = zeros(reso_num,1);
%         beta(index_b) = P\v;
%     end  
% end