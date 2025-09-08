function [R, stats] = montecarlo_moments_corr(t, lambda, Nsim)
% 使用蒙特卡洛方法估计归一化时间矩的相关系数矩阵
%
% 输入:
%   t      - K×1 时间向量
%   lambda - K×1 泊松强度向量
%   Nsim   - 模拟次数（建议 >= 10000）
%
% 输出:
%   R      - 3×3 相关系数矩阵，行列顺序为 [M1, M2, V]
%   stats  - 模拟出的 [M1, M2, V] 三列样本，Nsim × 3

% 矩阵准备
K = numel(t);
Lambda = sum(lambda);
stats = zeros(Nsim, 4);  % 每行存一组 [M1, M2, V]

for n = 1:Nsim
    % 对每个 bin 采样一个泊松数
    ni = poissrnd(lambda);   % K×1

    N_total = sum(ni);
    if N_total == 0
        % 跳过无探测事件的样本（避免除0）
        continue;
    end

    % 样本的时间统计量
    S1 = sum(t .* ni);        % ∑ t_i * n_i
    S2 = sum(t.^2 .* ni);     % ∑ t_i^2 * n_i

    M0 = N_total;
    M1 = S1 / N_total;              % 归一化一阶矩
    M2 = S2 / N_total;              % 归一化二阶矩
    V  = M2 - M1^2;                 % 方差

    stats(n, :) = [M0, M1, M2, V];
end

% 删除无效行（所有为 0 的跳过项）
stats(~any(stats,2), :) = [];

% 计算相关系数矩阵
R = corrcoef(stats);

end
