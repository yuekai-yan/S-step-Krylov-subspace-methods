function [Q, R] = GS(X, B)
    %Gram-Schmidt against other blocks computing R by backslash
    %Input:
    %X      n-by-m matrix, whose columns are orthonormal       
    %B      small block matrix of size n-by-s to be 
    %       orthogonalized with X
    %Output:
    %Q      orthogonalized matrix against X of size n-by-s
    %R      m-by-s matrix

    n = size(X, 1);
    m = size(X, 2);
    s = size(B, 2);
    Q = zeros(n, s);
    R = zeros(m, s);

    R = X \ B;      %m-by-s matrix
    Q = B - X * R;
end