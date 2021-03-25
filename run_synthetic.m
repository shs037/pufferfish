%% Compute the noise scale for MQMApprox, MQMExact, GK16 on synthetic data

epsilon = 1; % epsilon to achieve in pufferfish

% q = Pr(X_i = 0 | X_{i-1} = 0), p = Pr(X_i = 1 | X_{i-1} = 1)
% qs and ps initialize all possible q and p values
qs = 0.1:0.1:0.9;
ps = qs;

% initalize variables to save the noise parameters for each algorithm:
%   noise_mqm:        MQMExact
%   noise_mqm_orig:   MQMApprox
%   noise_inf:        GK16 (inferential privacy algorithm)
noise_mqm_exact = zeros(length(qs), length(ps)); noise_inf = zeros(length(qs), length(ps)); noise_mqm_approx = zeros(length(qs), length(ps));
% initalize variables to save the eigen-gap and pi_min
pis = zeros(length(qs), length(ps)); eigaps = zeros(length(qs), length(ps));

for T = [100] % T is the length of the chain
    for i_q = 1:length(qs)
        q = qs(i_q);
        for i_p = 1:length(ps)
            p = ps(i_p);
            
            % Run GK16
            gamma = inferentialPrivacy_cal_gamma(p, q);
            noise_inf(i_q, i_p) = 1 / inferentialPrivacy2(gamma, T, epsilon, 'paper exact'); %inferentialPrivacy(gamma, T)/(epsilon);
            
            % Run MQMApprox
            M = [1-q, q; 1-p, p]; % transition matrix
            eigap = eig(M); eigap = sort(abs(eigap), 'descend'); eigap = 1 - eigap(2); % eigap is the eigen-gap
            pi = M^100; pi = pi(1,:); % pi is (approximately) the stationary distribution
            pis(i_q, i_p) = min(pi); % save the eigen-gap and pi_min
            eigaps(i_q, i_p) = eigap;
            noise_mqm_approx(i_q, i_p) = k_findBest_2dir(pis(i_q, i_p), eigaps(i_q, i_p), epsilon, T);
            
            % Run MQMExact
            [c0, c1] = computeProbability(T, M, pi);
            [downstream, upstream1, upstream2] = exactRatioBinaryGenerate(M, T, c0, c1);
            noise_mqm_exact(i_q, i_p) = k_findBest_2dir_exact(downstream, upstream1, upstream2, epsilon, T);
        end
    end
end
