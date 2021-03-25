%% MQMExact for binary state, searching for all node X_i

function [cur_max, tmp, noise_min, best_a, best_b] = k_findBest_2dir_exact(downstream, upstream1, upstream2, epsilon, T)

noise_min = -ones(1, T); best_a = nan(1,T); best_b = nan(1,T);

% the worst way of adding noise
noise_largest = T / epsilon;

cur_max = -1;
for i = 2:T-1
    as = 1:i-1; % a as row
    bs = [1:T-i]'; % b as col
    
    delta_max = [downstream(bs, 1); 1]*[upstream1(as, 1)*upstream2(i,1); 1]';
    delta_min = [downstream(bs, 2); 1]*[upstream1(as, 2)*upstream2(i,2); 1]';
    
    delta = max(delta_max, 1./delta_min);
    delta = log(delta);
    
    tmp(i) = delta(1);
    
    lambda = epsilon - delta; lambda(lambda < 0) = nan; % lambda can't < 0
    
    noises = (repmat([as, i], length(bs)+1, 1) + repmat([bs; T-i+1], 1, length(as)+1) - 1) ./ lambda;
    
    [noise_min(i), best_idx] = min(noises(:));
    [best_b(i), best_a(i)] = ind2sub(size(noises), best_idx);
    if noise_min(i) > cur_max
        cur_max = noise_min(i);
    end
    
    if isnan(noise_min(i))
        cur_max = noise_largest;
        return;
    end
    
    if best_b(i) == T-i && best_a(i) == i-1
        disp('max_len is too small?')
    end
end

