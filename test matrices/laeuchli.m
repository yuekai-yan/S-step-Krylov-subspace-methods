function L = laeuchli(n, mu)
% LAEUCHLI Generate a Läuchli test matrix.
%   L = laeuchli(n, mu)              

    L = [ones(1,n); mu*eye(n-1, n)];
end
