function [Q, R] = CGS2(X, B)
    %Classic Gram-Schmidt with reorthogonalization against other blocks
    %Input:
    %X      n-by-m matrix, whose columns are Theta
    %       orthonormal
    %B      small block matrix of size n-by-s to be 
    %       Theta orthogonalized with X
    %Output:
    %Q      Theta orthogonalized matrix against X of size n-by-s
    %R      m-by-s matrix

    n = size(X, 1);
    m = size(X, 2);
    s = size(B, 2);
    Q = zeros(n, s);
    R = zeros(m, s);

    R0 = X' * B;  %m-by-s matrix
    Q0 = B - X * R0;
    %reorthogonalization
    R1 = X' * Q0;
    R = R0 + R1;    
    Q = Q0 - X * R1;
end