function [coef, res, SSR, yhat, sdr, ccoef, bcmean] = expregress2(x,y,sy,a,assertzero)
% EXPREGRESS2 is designed to optimize an exponential regression by
% minimizing the chi-squared of the residuals. 
%   INPUTS: 
%       x:  independent variables of data (vector)
%       y:  dependent variables of data, or measurements (vector)
%       sy: standard deviation of replicate measurements on y (vector)
%       a: parameter guesses. A vector where the a(1) is the
%       pre-exponential factor and a(2) is the rate coefficient
%       assertzero: binary flag where if==
%           1: use a(1) as a fixed parameter for the regression
%           0: use the a(1) as a true guess and let it be fit
%   OUTPUTS:
%       coef:   optimized coefficients, in the same order as a
%       res:    residuals of the regression model vs. data
%       r2:     R-squared of the model, or explained variance%
%       yhat:   estimated values of y using the model coefficients
%       sdr:    standard deviation of the residuals
%       ccoef:  Uncertainties in coefficient estimates
%           2x2 matrix [lower a(1), upper a(1); lower a(2), upper b(2)]
%       bcmean: Coefficients determined by bootstrapping 10^5 times.

% If the vectors are not columns, we will make them into columns.
if size(x, 2) ~= 1
    x = x';
end
if size(y,2) ~= 1
    y = y';
end
if size(sy,2) ~= 1
    sy = sy';
end

% Remove NaNs from the data
nans = isnan(x) | isnan(y);
x(nans) = [];
y(nans) = [];
sy(nans) = [];

% Default values for the coefficients
if nargin < 4
    a = [1, 0.01];
end

% Objective function to minimize
if assertzero
    k0 = a(2);
    C0 = a(1);
    chisq = @(k) sum( (y - C0.*exp(k.*x)).^2 ./ sy.^2);
    % Optimization
    k = fminsearch(chisq, k0);
    coef = [C0, k];
else
    chisq = @(a) sum( (y - a(1).*exp(a(2).*x)).^2 ./ sy.^2);
    % Optimization
    coef = fminsearch(chisq, a);
end



% Function to evaluate yhat
ya = @(x) coef(1).* exp(coef(2).*x);

yhat = ya(x);

% residuals
res = y-yhat;

% Goodness-of-fit metrics
n = size(x,1);
v = var(y);
UXSS = sum(res.^2);
SSR = UXSS;
SS = n*v;
XSS = SS-UXSS;
r2 = XSS/SS;
sdr = sqrt(XSS/n);

% Uncertainty in coefficient estimates. 
% This one is weird, so bear with me. We have vectors of y and sy, and from
% this I am going to make random numbers drawn from the distribution
% N(y(i),sy(i)). To refit the model, I need more uncertainties, but I am
% going to penalize the random numbers by adding sy(i) to the distance
% between the original mean (y(i)) and the estimate. This will give me a
% bootstrap set of parameter estimates. Their distribution will provide a
% 95% confidence interval. 

nb = 10000; % of bootstraps
rng('default')
pd = @(mu, sig) random(makedist("Normal", "mu",mu, "sigma",sig), [1,nb]);
yb = zeros(length(y), nb);
for ii = 1:size(y,1)
    yb(ii,:) = pd(y(ii),sy(ii));
end
syb = repmat(sy,1,nb)+(abs(yb-y));

ab = zeros(2, nb);
if assertzero
    ab(1,:) = C0;
    k0 = coef(2);
    C0 = coef(1);
    for ii = 1:nb
        chisq = @(k0) sum( (yb(:,ii) - C0.*exp(k0.*x)).^2 ./ syb(:,ii).^2);
        ab(2,ii) = fminsearch(chisq, k0)';
    end
else
    for ii = 1:nb
        chisq = @(a) sum( (yb(:,ii) - a(1).*exp(a(2).*x)).^2 ./ syb(:,ii).^2);
        ab(:,ii) = fminsearch(chisq, coef)';
    end
end

% So, after forcing MATLAB to do all this estimation, might as well output
% the parameters that got bootstrapped to hell and back.
if sum(sum(isnan(ab))) == size(ab,1)*size(ab,2)
    bcmean = [NaN,NaN];
    ccoef = nans(2,2);
    disp("Fitting failed, imaginary or NaN parameters. Setting to NaN.")
else
    bcmean = [mean(ab(1,:)),-exp(mean(log(-ab(2,:))))];

    ccoef = zeros(2,2);
    ccoef(1,:) = norminv([0.05, 0.95], bcmean(1), std(ab(1,:)));
    ccoef(2,:) = -logninv([0.05, 0.95], mean(log(-ab(2,:))), std(log(-ab(2,:))));
end

% figure
% subplot(3,1,1)
% histogram(ab(1,:))
% hold on
% plot(ccoef(1,:),[1,1]*500)
% subplot(3,1,2)
% histogram(ab(2,:))
% hold on
% plot(ccoef(2,:),[1,1]*500)
% subplot(3,1,3)
% histogram(log(-ab(2,:)))
% Those figures were to show that the parameter a(1) is normally
% distributed and a(2) is lognormal. Hence, we calculate the CI a little
% differently. 
end