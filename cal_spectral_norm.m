function res = cal_spectral_norm(A)

res = sqrt(max(eig(A'*A)));