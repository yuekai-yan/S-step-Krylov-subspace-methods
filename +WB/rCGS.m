function [Q, R] = rCGS(X, Theta)
    %Randomized Classic Gram-Schmidt within a block
    %Input: 
    %X n-by-m matrix to be Theta orthogonalized
    %Theta  random sketch of size d-by-n
    %Output: 
    %Q n-by-m matrix, column-oriented Theta-orthonormal matrix
    %R m-by-m matrix, coefficients of the basis spanned by Q

    n = size(X, 1);
    m = size(X, 2);
    d = size(Theta, 1);
    Q = zeros(n, m);
    R = zeros(m, m);
    S = zeros(d, m);
    P = Theta * X;  %d-by-m matrix
    for i = 1:m
        w = S(:, 1:i-1)' * P(:, i); %(i-1)-by-1 matrix
        R(1:i-1, i) = w;
        q = X(:, i) - Q(:, 1:i-1) * w;
        s = Theta * q;
        R(i, i) = norm(s);
        S(:, i) = s / R(i, i);
        Q(:, i) = q / R(i, i);
    end
end