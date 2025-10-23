function [Q, R] = rMGS(X, Theta)
    %Randomized Modified Gram-Schmidt within a block
    %Input:
    %X n-by-m matrix to be orthogonalized
    %Theta  random sketch of size d-by-n
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
        q = X(:, i);
        p = P(:, i);
        for j = 1:i-1
            w = S(:, j)' * p;
            R(j, i) = w;
            q = q - Q(:, j) * w;
        end
        s = Theta * q;
        R(i, i) = norm(s);
        S(:, i) = s / R(i, i);
        Q(:, i) = q / R(i, i);
    end
end