%% Run GK16 algorithm
% Input:
%   gamma: gamma value in the algorithm
%   T: length of the Markov chain
%   nu_pfp: the pufferfish privacy parameter we aim to achieve
%   option: which algorithm to use

function dp_epsilon = inferentialPrivacy2(gamma, T, nu_pfp, option)
Gamma = diag(gamma*ones(1, T-1), -1) + diag(gamma*ones(1, T-1), 1);

if gamma > 1/2
    if cal_spectral_norm(Gamma) >= 1
            disp('spectral norm larger than 1')
        dp_epsilon = nan;
        return;
    end
end

if nargin == 3
    option = 'paper approx';
end

switch option
    case 'paper exact'
        Phi = inv(eye(T) - Gamma);
        dp_epsilon = nu_pfp / max(2*sum(Phi,2));
    case 'paper approx'
        dp_epsilon = nu_pfp*(1-2*gamma*T) / 2;
    case 'tight'
        dp_epsilon = max(nu_pfp*(1-2*gamma*T)/(2-2*gamma*T), (nu_pfp - 2*gamma*T)/2);        
end
dp_epsilon = max(0, dp_epsilon);
