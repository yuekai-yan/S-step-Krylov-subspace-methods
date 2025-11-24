function [Q, R] = rMGS(X, B, Theta)
    %Randomized modified Gram-Schmidt against other blocks
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
    Q = B;
    R = zeros(m, s);

    P = Theta * B;  %d-by-s matrix
    S = Theta * X;  %d-by-m matrix
    for i = 1:m
        W = S(:, i)' * P;  %1-by-s vector
        R(i, :) = W;
        P = P - S(:, i) * W;
        Q = Q - X(:, i) * W;
    end
end