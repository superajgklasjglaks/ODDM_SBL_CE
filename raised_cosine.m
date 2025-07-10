% function y = raised_cosine(t, beta, A, T)
%     % raised_cosine - 计算升余弦函数在给定点的值
%     %
%     % 输入参数:
%     % t - 时间点
%     % beta - 滚降系数 (0 <= beta <= 1)
%     % A - 幅度（函数的峰值）a
%     % T - 符号间隔（函数的宽度，通常为符号周期）
%     %
%     % 输出参数:
%     % y - 升余弦函数在 t 点的值
% 
%     % 检查滚降系数是否在有效范围内
%     if beta < 0 || beta > 1
%         error('滚降系数 beta 必须在 0 和 1 之间');
%     end
% 
%     % 检查 T 是否为正数
%     if T <= 0
%         error('符号间隔 T 必须为正数');
%     end
% 
%     % 归一化时间
%     t_norm = t / T;
% 
%     % 计算升余弦函数
%     if t == 0
%         y = A;
%     elseif abs(t_norm) == 1 / (2 * beta)
%         y = A * (sin(pi / (2 * beta))) / (pi / (2 * beta));
%     else
%         y = A * (sin(pi * t_norm) / (pi * t_norm)) * (cos(beta * pi * t_norm) / (1 - (2 * beta * t_norm)^2));
%     end
% end

function h = raised_cosine(t, T, alpha)
% 升余弦滤波器时域函数
% 输入:
%   t: 时间向量
%   T: 符号周期
%   alpha: 滚降因子 (0 ≤ alpha ≤ 1)
% 输出:
%   h: 升余弦滤波器的时域响应

    h = zeros(size(t));
    idx = (t ~= 0); % 排除 t=0 的点
    
    % 主公式计算
    numerator = sin(pi * t(idx) / T) .* cos(pi * alpha * t(idx) / T);
    denominator = (pi * t(idx) / T) .* (1 - (2 * alpha * t(idx) / T).^2);
    h(idx) = numerator ./ denominator;
    
    % 处理 t=0 的极限 (h(0) = 1)
    h(t == 0) = 1;
    
    % 处理 t = ±T/(2alpha) 的极限
    t_special_pos = T / (2 * alpha);
    t_special_neg = -T / (2 * alpha);
    if alpha ~= 0
        h(abs(t - t_special_pos) < 1e-10) = (alpha/2) * sin(pi/(2*alpha));
        h(abs(t - t_special_neg) < 1e-10) = (alpha/2) * sin(pi/(2*alpha));
    end
end