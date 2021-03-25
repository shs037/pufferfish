%% MQMExact for multi-state, searching for all node X_i

function [cur_max, tmp, noise_min, best_a, best_b] = k_findBest_2dir_exact_multi(downstream, upstream1, upstream2, epsilon, T)

noise_min = -ones(1, T); best_a = nan(1,T); best_b = nan(1,T);

% the worst way of adding noise
noise_largest = T / epsilon;

cur_max = -1;
for i = 2:T-1
% for i = T/2
    as = 1:i-1; % a as row
    bs = [1:T-i]'; % b as col
    
    delta = -ones(length(bs)+1, length(as)+1);
    for j = 1:size(downstream, 2)
        delta = max(delta, [downstream(bs, j); 1]*[upstream1(as, j)*upstream2(i, j); 1]');
    end
    
    delta = log(delta);
    
    tmp(i) = delta(1);
    
    lambda = epsilon - delta; lambda(lambda < 0) = nan; % lambda can't < 0
    
    noises = (repmat([as, i], length(bs)+1, 1) + repmat([bs; T-i+1], 1, length(as)+1) - 1) ./ lambda;
    
    [noise_min(i), best_idx] = min(noises(:));
    [best_b(i), best_a(i)] = ind2sub(size(noises), best_idx);
    %     disp([best_b(i), best_a(i)])
    if noise_min(i) > cur_max
        cur_max = noise_min(i);
    end
    
    if isnan(noise_min(i))
        cur_max = noise_largest;
        return;
        %         error(['at i=' num2str(i) ', cannot find valid a/b. epsilon should be larger']);
    end
    
    if best_b(i) == T-i && best_a(i) == i-1
        disp('max_len is too small?')
    end
end

