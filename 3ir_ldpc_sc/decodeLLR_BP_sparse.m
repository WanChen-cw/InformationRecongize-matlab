function [vHat, n] = decodeLLR_BP_sparse(llr, H, iteration, s)
% 获取校验矩阵H的大小
[M, N] = size(H);

% 先验概率对数似然值
Lci = llr;

% E校验矩阵Lr
Lrji = sparse(M, N);

% M变量矩阵初始化
Lqij = sparse(H) .* repmat(Lci, M, 1);

n = 0;
while n < iteration
    % 译码轮次
    fprintf('Iteration: %d\n', n);
    
    % 变量节点矩阵的符号信息
    alphaij = sign(Lqij);
    % 变量节点矩阵的绝对值信息
    betaij = abs(Lqij);
    
    % 校验节点矩阵的更新
    idx = find(H);
    [i, c1] = ind2sub([M, N], idx);
    Lrji(idx) = sparse(i, c1, Lrji(idx), M, N);
    Lrji(idx) = Lrji(idx) + sparse(i, c1, betaij(i, c1), M, N);
    Lrji(idx) = Lrji(idx) - sparse(i, c1, alphaij(i, c1) .* betaij(i, c1), M, N);
    for k = 1:N
        if s(k) == 1
            Lrji(:, k) = -Lrji(:, k);
        end
    end
    
    % 变量节点矩阵的更新
    Lqij = sparse(H) .* repmat(Lci, M, 1) - Lrji;
    
    % 计算后验概率
    LQ = sum(Lrji, 1) + Lci;
    
    % 进行硬判决
    vHat = (LQ < 0);
    
    % 判断译码是否成功
    if isequal(mod(vHat * H', 2), s)
        break;
    else
        n = n + 1;
    end
end % while
