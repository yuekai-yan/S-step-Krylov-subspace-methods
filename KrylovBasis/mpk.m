function [V, B] = mpk(A, q, s, beta, gamma)
    %Inputs: 
    %  A          matrix of size (n, n)
    %  q          start vector of size (n, 1)
    %  s          the dimension of Krylov subspace
    %  beta       vector of size (s, 1)
    %  gamma      vector of size (s, 1)
    %
    %Outputs:
    %  V          basis of size (n, s)
    %  B          change-of-basis matrix such that 
    %             A * V(:, 1:s) = V(:, 1:(s+1)) * B
    
    if nargin < 4
        beta = zeros(s, 1);
    end

    if nargin < 5
        gamma = ones(s, 1);
    end
    n = size(A, 1);
    V_1 = zeros(n, s+1); %The matrix is augmented with an additional
                         %column p_s(A) * q
    B = zeros(s+1, s);
    P = eye(n);

    v = q / norm(q);
    V_1(:, 1) = v;

    for j = 1:s
        P = (A - beta(j) * eye(n)) * P;
        V_1(:, j+1) = P * v;
    end

    V = V_1(:, 1:s);
    B = V_1 \ (A * V);
end