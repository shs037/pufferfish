function [downstream, upstream1, upstream2] = exactRatioBinaryGenerate(M, T, c0, c1)
% downstream: length T-1 vector. downstream(b): max_{x_{i+b}} P(x_{i+b} | X_i = 0)/P(x_{i+b} | X_i = 1) (this is the same for all i)
% upstream_1: length T-1 vector. upstream_1(a): max_{x_{i-a}} P(X_i = 0 | x_{i-a})/P(X_i = 1 | x_{i-a}) (this is the same for all i)
% upstream_2: length T vector.   upstream_2(i): P(X_i=1)/P(X_i=0) (different for different i)

k = size(M, 1);
if k ~= 2
    error('this handles 2 by 2 transition matrix only')
end

Mpower = zeros(k, k, T-1); Mpower(:,:,1) = M;
for i = 2:T-1
    Mpower(:,:,i) = M*Mpower(:,:,i-1);
end

% downstream
downstream = [permute(max(Mpower(1,:,:)./Mpower(2,:,:), [], 2), [3, 1, 2]), ...
              permute(min(Mpower(1,:,:)./Mpower(2,:,:), [], 2), [3, 1, 2])];

% upstream
upstream1 = [permute(max(Mpower(:,1,:)./Mpower(:,2,:), [], 1), [3, 1, 2]), ...
             permute(min(Mpower(:,1,:)./Mpower(:,2,:), [], 1), [3, 1, 2])];

if nargin == 4 % given fixed c0, c1
    upstream2 = repmat(c1./c0, 1, 2);
elseif nargin == 2 % take max c0, c1
    upstream2 = [1./upstream1(:,2), 1./upstream1(:,1)];
    upstream2 = [nan, nan; upstream2];
else
    error('exactRatioBinaryGenerate input formate error')
end
    
