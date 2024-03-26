function chisq = chisq_quantum_invariant(mc,sc, I,eps,alph,Q,t,lambda,r)
% Gives a chi-squared value for the values from eval_left.m (data) and
% int_quantum.m (estimate)

kp = -log(mc(2:end)./mc(1))./t(2:end);
bad = find(isnan(kp) | isinf(kp) | kp<0);
kp(bad, :) = []; t(bad, :) = [];
mc(bad+1,:) = []; sc(bad + 1, :) = [];
I( :, bad+1) = [];

[y, sy] = eval_left(mc,sc);

Y = 1e9*int_quantum(mc, I',eps,alph,Q,t,lambda,r);

chisq = sum((y - Y').^2 ./ (sy.^2));