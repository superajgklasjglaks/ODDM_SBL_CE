function h_prime = raised_cosine_derivative(t, T, alpha)
% 升余弦滤波器的导数
% 输入:
%   t: 时间向量
%   T: 符号周期
%   alpha: 滚降因子 (0 ≤ alpha ≤ 1)
% 输出:
%   h_prime: 升余弦滤波器的时域导数

    h_prime = zeros(size(t));
    idx = (t ~= 0) & (abs(t) ~= T/(2*alpha)); % 排除奇点
    
    % 符号导数公式的直接实现
    term1 = (pi * t(idx)) / T .* cos(pi * t(idx) / T) - sin(pi * t(idx) / T);
    term1 = term1 ./ (pi * t(idx) / T .* sin(pi * t(idx) / T));
    
    term2 = -alpha * tan(pi * alpha * t(idx) / T);
    
    term3 = (8 * alpha^2 * t(idx)) / (pi * T);
    term3 = term3 ./ (1 - (2 * alpha * t(idx) / T).^2);
    
    % 组合三项
    h_t = raised_cosine(t(idx), T, alpha);
    h_prime(idx) = h_t .* (term1 + term2 + term3);
    
    % 处理奇点 (t=0 和 t=±T/(2alpha) 的导数需单独定义)
    % (1) t=0 的导数: 通过对称性可知 h'(0)=0
    h_prime(t == 0) = 0;
    
    % (2) t=±T/(2alpha) 的导数: 需用极限计算
    % if alpha ~= 0
    %     t_special = T / (2 * alpha);
    %     [~, idx_pos] = min(abs(t - t_special));
    %     [~, idx_neg] = min(abs(t + t_special));
    % 
    %     % 使用数值差分逼近极限
    %     dt = 1e-6;
    %     h_prime(idx_pos) = (raised_cosine(t_special + dt, T, alpha) - ...
    %                        raised_cosine(t_special - dt, T, alpha)) / (2*dt);
    %     h_prime(idx_neg) = (raised_cosine(-t_special + dt, T, alpha) - ...
    %                        raised_cosine(-t_special - dt, T, alpha)) / (2*dt);
    % end
end