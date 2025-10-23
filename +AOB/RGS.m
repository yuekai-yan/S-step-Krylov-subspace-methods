function [Q, R] = RGS(X, B, Theta)
    %Randomized Gram-Schmidt against other blocks
    %Input:
    %X      n-by-m matrix, whose columns are Theta
    %       orthonormal 
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
    S = Theta * X;  %d-by-m matrix
    R = S \ P;      %m-by-s matrix
    Q = B - X * R;
end