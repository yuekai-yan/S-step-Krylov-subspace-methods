function [Q, R] = Whitening(X)
    %Whitening within a block
    %Input: 
    %X n-by-m matrix to be Theta orthogonalized
    %Output: 
    %Q n-by-m matrix, column-oriented Theta orthonormal matrix
    %R m-by-m matrix, coefficients of the basis spanned by Q

    n = size(X, 1);
    m = size(X, 2);
    Q = zeros(n, m);
    R = zeros(m, m);
    [Q, R] = qr(X, 0);
end