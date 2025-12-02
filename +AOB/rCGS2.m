function [Q, R] = rCGS2(X, S, B, Theta)
    %Randomized classic Gram-Schmidt with reorthogonalization 
    %against other blocks
    %Input:
    %X      n-by-m matrix, whose columns are Theta
    %       orthonormal
    %S      d-by-m matrix, sketched matrix of X
    %B      small block matrix of size n-by-s to be 
    %       Theta orthogonalized with X
    %Theta  d-by-n matrix
    %Output:
    %Q      Theta orthogonalized matrix against X of size n-by-s
    %R      m-by-s matrix

    n = size(X, 1);
    m = size(X, 2);
    s = size(B, 2);
    Q = zeros(n, s);
    R = zeros(m, s);

    P = Theta * B;  %d-by-s matrix
    %S = Theta * X;  %d-by-m matrix
    R0 = S' * P;  %m-by-s matrix
    Q0 = B - X * R0;
    %reorthogonalization
    S0 = Theta * Q0;
    R1 = S' * S0;
    R = R0 + R1;
    Q = Q0 - X * R1;
end