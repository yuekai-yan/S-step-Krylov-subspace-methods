function [Q, R] = rMGS(X, S, B, Theta)
    %Randomized modified Gram-Schmidt against other blocks
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
    R = zeros(m, s);

    P = Theta * B;  %d-by-s matrix
    %S = Theta * X;  %d-by-m matrix
    for i = 1:m
        W = S(:, i)' * P;  %1-by-s vector
        R(i, :) = W;
        P = P - S(:, i) * W;
    end
    Q = B - X * R;
end