function [Q, R] = MGS(X)
    %Modified Gram-Schmidt within a block
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
        q = X(:, i);
        for j = 1:i-1
            w = Q(:, j)' * q;
            R(j, i) = w;
            q = q - Q(:, j) * w;
        end
        R(i, i) = norm(q);
        Q(:, i) = q / R(i, i);
    end
end