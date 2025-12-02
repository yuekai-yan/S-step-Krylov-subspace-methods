function [Q, R] = GS(X)
    %Gram-Schmidt within a block computing R by backslash
    %Input: 
    %X        n-by-m matrix to be Theta-orthogonalized
    %Output: 
    %Q        n-by-m matrix, column-oriented orthonormal matrix
    %R        m-by-m matrix, upper triangular factor

    n = size(X, 1);
    m = size(X, 2);
    Q = zeros(n, m);
    %store the coefficients of the basis spanned by S
    R = zeros(m, m);
    for i = 1:m
        %solve the sketched least square problem
        R(1:i-1, i) = Q(:, 1:i-1) \ X(:, i);
        q = X(:, i) - Q(:, 1:i-1) * R(1:i-1, i);
        R(i, i) = norm(q);
        Q(:, i) = q / R(i, i);
    end
end