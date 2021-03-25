%% Compute the noise scale of MQMApprox, MQMExact and GroupDP on real data and plot the histogram.
% Input:
%   c:          original chain
%   M:          transition matrix

epsilon = 1; % privacy parameter for Pufferfish

%  Calculate eigen-gap (eigap) and pi_min (pi)
eigap = eig(M*reversibleChain_reverse(M)); eigap = sort(abs(eigap), 'descend'); eigap = (1 - eigap(2))/2;
pi = M^100; pi = pi(1,:); pi = min(pi);

% max range to search for Markov Quilt. Can also be set to fixed value.
l = 2 * 2 * 2*ceil(log((exp(epsilon/6) + 1) / (exp(epsilon/6) -1) / pi) / (2*eigap));

% Calculate noise scale for MQMApprox (need to add Lap(2*noise_mqm_approx) to each histogram bar)
noise_mqm_approx = k_findBest_2dir(pi, eigap, epsilon, l);

% Calculate noise scale for MQMExact (need to add Lap(2*noise_mqm_exact) to each histogram bar)
[downstream, upstream1, upstream2] = exactRatioMultiGenerate(M, l);
noise_mqm_exact = k_findBest_2dir_exact_multi(downstream, upstream1, upstream2, epsilon, l);

% Calculate noise scale for GroupDP (need to add Lap(2*noise_groupdp) to
% each histogram bar). If c contains several subchains, then length(c)
% should be replaced by the max-length of the subchains.
noise_groupdp = length(c) / epsilon;



% Plot the results
figure(1); hold all;
plotHist(c); % plot the original histogram
plotHistNoise(noise_groupdp, c, -.3, 4, '+', 20); % plot group-dp for 20 times
plotHistNoise(noise_mqm_approx, c, 0, 5, 'x', 20); % plot MQMApprox for 20 times 
plotHistNoise(noise_mqm_exact, c, .3, 6, 'o',20); % plot MQMExact for 20 times


