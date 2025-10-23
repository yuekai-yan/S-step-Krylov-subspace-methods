function [Q, R] = rCGS2(X, Theta)
    %Randomized Classic Gram-Schmidt with reorthogonalization within a
    %block
    %Input: 
    %X n-by-m matrix to be orthogonalized
    %Theta random sketch of size (d, n)
    %Output: 
    %Q n-by-m matrix, column-oriented orthonormal matrix
    %R m-by-m matrix, coefficients of the basis spanned by Q

    n = size(X, 1);
    m = size(X, 2);
    d = size(Theta, 1);
    Q = zeros(n, m);
    S = zeros(d, m);
    R = zeros(m, m);
    P = Theta * X;  %d-by-m matrix
    for i = 1:m
        w0 = S(:, 1:i-1)' * P(:, i);  %(i-1)-by-1 vector
        q0 = X(:, i) - Q(:, 1:i-1) * w0;
        s0 = Theta * q0;
        %reorthogonalization
        w = S(:, 1:i-1)' * s0;
        R(1:i-1, i) = w0 + w;
        q = q0 - Q(:, 1:i-1) * w;

        s = Theta * q;
        R(i, i) = norm(s);
        S(:, i) = s / R(i, i);
        Q(:, i) = q / R(i, i);
    end
end