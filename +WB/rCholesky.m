function [Q, R] = rCholesky(X, Theta)
    %Randomized Cholesky within a block
    %Input: 
    %X n-by-m matrix to be Theta orthogonalized
    %Theta  random sketch of size d-by-n
    %Output: 
    %Q n-by-m matrix, column-oriented Theta orthonormal matrix
    %R m-by-m matrix, coefficients of the basis spanned by Q

    n = size(X, 1);
    m = size(X, 2);
    d = size(Theta, 1);
    Q = zeros(n, m);
    R = zeros(m, m);
    S = zeros(d, m);
    P = Theta * X;  %d-by-m matrix
    [S, R] = qr(P, 0);
    Q = X / R;
end