function [aod_bar, aoa_bar] = MIMO_First_Order_Linear_Approximation(N_aod, M_aoa, aod_max, r_aod, r_aoa)
% Massive MIMO 一阶线性近似函数
% 模仿ODDM中的First_Order_Linear_Approximation
% 输入:
%   N_aod: 发射角网格数
%   M_aoa: 接收角网格数
%   aod_max: 发射角最大值
%   r_aod: 发射角分辨率
%   r_aoa: 接收角分辨率
% 输出:
%   aod_bar: 发射角网格
%   aoa_bar: 接收角网格

virtual_size = N_aod * M_aoa;
aod_bar = zeros(virtual_size, 1);
aoa_bar = zeros(virtual_size, 1);

% 生成发射角网格
aod_grid = linspace(-aod_max, aod_max, N_aod);
% 生成接收角网格
aoa_grid = linspace(-aod_max, aod_max, M_aoa);  % 假设接收角范围与发射角相同

% 填充网格
index = 1;
for n = 1:N_aod
    for m = 1:M_aoa
        aod_bar(index) = aod_grid(n);
        aoa_bar(index) = aoa_grid(m);
        index = index + 1;
    end
end

end