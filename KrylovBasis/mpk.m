function [V, B] = mpk(A, q, s, beta, gamma)
    %Inputs: 
    %A          function handle: (n×1 vector) -> (n×1 vector)
    %q          start vector of size (n, 1)
    %s          the dimension of Krylov subspace
    %beta       vector of size (s, 1)
    %gamma      vector of size (s, 1)
    %
    %Outputs:
    %V          basis of size (n, s)
    %B          change-of-basis matrix of size (s+1, s) such that 
    %             A * V_1(:, 1:s) = V_1(:, 1:(s+1)) * B
    
    if nargin < 4
        beta = zeros(s, 1);
    end

    if nargin < 5
        gamma = ones(s, 1);
    end

    n = length(q);
    V_1 = zeros(n, s+1);
    B = zeros(s+1, s);

    v = q / norm(q);
    V_1(:, 1) = v;

    for j = 1:s
        w = A(V_1(:, j)) - beta(j) * V_1(:, j);
        V_1(:, j+1) = w / gamma(j);
    end

    V = V_1(:, 2:s+1);
    
    B(1:s, 1:s) = diag(beta);
    B(2:(s+1), 1:s) = B(2:(s+1), 1:s) + diag(gamma);
end