function [Q, R] = RGS(X, Theta)
    %Randomized Gram-Schmidt within a block
    %Input: 
    %X        n-by-m matrix to be Theta-orthogonalized
    %Theta    random sketch of size d-by-n
    %Output: 
    %Q        n-by-m matrix, column-oriented Theta-orthonormal matrix
    %R        m-by-m matrix, upper triangular factor

    n = size(X, 1);
    m = size(X, 2);
    d = size(Theta, 1);
    Q = zeros(n, m);
    %store the coefficients of the basis spanned by S
    R = zeros(m, m);
    %store the sketched orthonormal basis
    S = zeros(d, m);
    for i = 1:m
        p = Theta * X(:, i);
        %solve the sketched least square problem
        R(1:i-1, i) = S(:, 1:i-1) \ p;
        q = X(:, i) - Q(:, 1:i-1) * R(1:i-1, i);
        s0 = Theta * q;
        R(i, i) = norm(s0);
        S(:, i) = s0 / R(i, i);
        Q(:, i) = q / R(i, i);
    end
end