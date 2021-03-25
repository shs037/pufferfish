%% Calculate gamma for GK16 algorithm
function [gamma,f1,f2,f3,f4] = inferentialPrivacy_cal_gamma(p, q)

f1 = (p^2 + p*(1-q))/(p^2 + q*(1-p));
f2 = (p/q)*((1-q)^2+q*(1-p)) / ((1-p)*(1-q) + p*(1-p));
f3 = (1-p)/(1-q) * ((1-q)*q + q*p) / ((1-p)*q+p^2);
f4 = ((1-q)^2 + p*(1-q))/((1-q)^2 + q*(1-p));

tmp = [f1, f2, f3, f4];
gamma = max([tmp, 1./tmp]);

gamma = log(gamma)/2;