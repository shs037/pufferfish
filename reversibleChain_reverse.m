%% Generate the reverse of matrix p
function [pp] = reversibleChain_reverse(p)

nState = length(p);

pi = p^4096; pi = pi(1,:)';
pp = repmat(pi, 1, nState) .* p ./ repmat(pi', nState, 1); pp = pp';

end