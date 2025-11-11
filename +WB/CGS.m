function [Q, R] = CGS(X)
    %Classic Gram-Schmidt within a block
    %Input: 
    %X n-by-m matrix to be orthogonalized
    %Output: 
    %Q n-by-m matrix, column-oriented orthonormal matrix
    %R m-by-m matrix, coefficients of the basis spanned by Q

    n = size(X, 1);
    m = size(X, 2);
    Q = zeros(n, m);
    R = zeros(m, m);

    for i = 1:m
        w = Q(:, 1:i-1)' * X(:, i); %(i-1)-by-1 matrix
        R(1:i-1, i) = w;
        q = X(:, i) - Q(:, 1:i-1) * w;
        R(i, i) = norm(q);
        Q(:, i) = q / R(i, i);
    end
end