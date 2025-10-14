function RitzValues = getRitzValues(Afun, x, s, reorth_tol)
% Input:
%   Afun       --  function handle that effects a matrix-vector multiplication,
%                   e.g., Afun = @(x) A*x if A is a matrix.
%   x          --  starting vector
%   s          --  step size
%   reorth_tol --  Reorthogonalization tolerance (default 0.7).  A
%                  tolerance of 0.0 implies *no* reorthogonalization.
%
% Output:
%   RitzValues --  s Ritz values of A

if nargin < 4 || isempty(reorth_tol)
    reorth_tol = 0.7;
end

u = x / norm(x);
U = zeros(length(x), s+1);
U(:, 1) = u;
H0 = zeros(s+1, s);
for k = 1:s
    w = Afun(u);
    H0(1:k, k) = U(:,1:k)' * w;
    u = w - U(:,1:k) * H0(1:k, k);
    if norm(u) < reorth_tol * norm(w)
        h_corr = U(:,1:k)' * u;
        H0(1:k, k) = H0(1:k, k) + h_corr;
        u = u - U(:,1:k) * h_corr; 
    end
    H0(k+1, k) = norm(u);
    u = u / H0(k+1, k);
    U(:,k+1) = u;
end
H = H0(1:s, 1:s);
RitzValues = eig(H);
end