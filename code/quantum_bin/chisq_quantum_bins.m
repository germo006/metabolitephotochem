function chisq = chisq_quantum_bins(mc,sc, W,eps,Qlength, Qbin,Qbounds, t,lambda)
% Gives a chi-squared value for the values from eval_left.m (data) and
% int_quantum.m (estimate)

Q = ones(Qlength,1);
for ii=1:size(Qbounds,1)-1
    Q(Qbounds(ii):Qbounds(ii+1)) = Qbin(ii);
end

t(1,:) = [];
kp = -log(mc(2:end)./mc(1))./t;
bad = find(isnan(kp) | isinf(kp) | kp<0);
kp(bad, :) = []; t(bad, :) = [];
mc(bad+1,:) = []; sc(bad + 1, :) = [];
W(bad, :) = [];

[y, sy] = eval_left(mc,sc);
Y = zeros(size(kp));
for ii=1:length(kp)
    Y(ii) = int_quantum(W(ii,:)',eps,Q,kp(ii),t(ii),lambda);
end

chisq = sum((y - Y).^2 ./ (sy.^2));