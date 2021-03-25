%% MQMApprox, searching from the middle node
% Input:
%   pi_min, eig_gap: pi_min and eigen-gap of the transition matrix
%   epsilon: pufferfish privacy parameter
%   T: length of the chain
% 
% Return: (only the first variable is actually used in the algorithm)
%   cur_max: the final noise scale
%   cur_len: length of the final Markov Quilt
%   cur_ab: (a, b) pair of the final Markov Quilt
%   noise_min: noise scale for each X_i
%   best_a, best_b: best (a, b) for each X_i
%   i: stopped search at this X_i

function [cur_max, cur_len, cur_ab, noise_min, best_a, best_b, i] = k_findBest_2dir(pi_min, eig_gap, epsilon, T)
if isnan(epsilon)
    cur_max = nan;
    return;
end

noise_min = -nan(1, T-1); best_a = noise_min; best_b = noise_min;

cur_max = -1; cur_len = -1; cur_ab = [-1 -1]; cur_status = '';
i = floor(T / 2);
while(true)
    as = 1:i-1; % a as row
    bs = [1:T-i]'; % b as col
    
    Delta_a = 1/pi_min * exp(-eig_gap*as); Delta_a(Delta_a >= 1) = nan; % nan if denominator (1-Delta) is < 0
    Delta_b = 1/pi_min * exp(-eig_gap*bs); Delta_b(Delta_b >= 1) = nan;
    
    delta_a = log((1+Delta_a)./(1-Delta_a)); % ratio automatically > 1
    delta_b = log((1+Delta_b)./(1-Delta_b));
    delta_a = [delta_a, 0]; delta_b = [delta_b; 0]; % handle the boundary (the if in FBB)
    delta = 2*repmat(delta_a, length(bs)+1, 1) + repmat(delta_b, 1, length(as)+1);
    
    lambda = epsilon - delta; lambda(lambda < 0) = nan; % lambda can't < 0
    
    noises = (repmat([as, i], length(bs)+1, 1) + repmat([bs; T-i+1], 1, length(as)+1) - 1) ./ lambda;
    
    [noise_min(i), best_idx] = min(noises(:));
    [best_b(i), best_a(i)] = ind2sub(size(noises), best_idx);
    if noise_min(i) > cur_max
        cur_max = noise_min(i);
        cur_len = best_a(i)+best_b(i)-1;
        cur_ab = [best_a(i), best_b(i)];
        if abs(cur_len / lambda(best_idx) - cur_max) > eps
            disp('error')
        end
        
    end
    
    if isnan(noise_min(i))
        error('cannot find valid a/b. epsilon should be larger');
    end
    
    if best_b(i) == T-i+1 && best_a(i) == i
        break;
    elseif best_b(i) <= T - i && best_a(i) <= i-1
        break;
    elseif best_a(i) == i % no left
        if strcmp(cur_status, 'no right')
            break;
        end
        cur_status = 'no left';
        i = i+1;
    elseif best_b(i) == T - i + 1
        if strcmp(cur_status, 'no left')
            break;
        end
        cur_status = 'no right';
        i = i-1;
    else
        disp([best_a(i), best_b(i), i, T])
        error('if-else problem.')
    end
%     disp('new round')
end

% the worst way of adding noise
noise_largest = T / epsilon;
if abs(cur_max - noise_largest) < eps
%     disp('same as global sensitivity * T');
end

% plot(2:max_len-1, noise_min(2:end))
% figure(2)
% hold all
% plot(2:max_len-1, best_b(2:end) - best_a(2:end))
