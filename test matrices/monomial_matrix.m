function M = monomial_matrix(x, d)
if nargin < 2 
    d = numel(x)-1; 
end
x = x(:);
M = x .^ (0:d);
end