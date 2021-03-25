%% Compute Pr(X_i = 1) and Pr(X_i = 0) for each X_i given initial distribution ini

function [c0, c1] = computeProbability(T, M, ini)

k = size(M, 1);
c = zeros(T, k);
c(1, :) = ini;
for i = 2:T
    c(i, :) = c(i-1, :)*M;
end
c0 = c(:, 1);
c1 = c(:, 2);
