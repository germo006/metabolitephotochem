function [y, sy] = eval_left(mc,sc)
% This function takes successive measurements of tryptophan (let's not kid
% ourselves like I'm doing anything else) and calculates the left-hand side
% of the function I need to minimize. It's the "data" instead of the
% estimate. 
%
%   INPUTS
%       mc: nx1 mean concentrations, nM
%       sc: nx1 standard deviations of measurements, nM
%   
%   OUTPUTS
%       y:  n-1x1 vector of C0-Ct / 2.3C0 for each non-initial point
%       sy: n-1x1 vector of errors on this calculation 
%

y = (mc(1) - mc(2:end)) ./ (2.303*mc(1));
sy = y.*sqrt((sc(1)/mc(1)).^2 + (sc(2:end)./mc(2:end)).^2);