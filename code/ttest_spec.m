function h = ttest_spec(x,y, alph)
%TTEST_SPEC is meant to run down two triplicate sets of measurements and
%perform a two-tailed t-test at each row. 
%   INPUT: 
%       x, an n x 3 matrix
%       y, an n x 3 matrix
%       alph, significance level. default is 0.01
%   OUTPUT
%       h, an n x 1 vector of 0s and 1s, where 1 is a place where the
%       measurements in row i of x differ from y at a significance of alph
%

if(~exist("alph", "var"))
    alph = 0.01;
end

h = zeros(size(x,1),1);
for ii = 1:size(x,1)
    h(ii) = ttest2(x(ii,:), y(ii,:), "Alpha", alph, "Tail", "both", "Vartype", "equal");
end 

h(isnan(h)) = 0;
h = logical(h);

end

