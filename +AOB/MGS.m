function [Q, R] = MGS(X, B)
    %Randomized modified Gram-Schmidt against other blocks
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
    Q = B;

    for i = 1:m
        W = X(:, i)' * B;  %1-by-s vector
        R(i, :) = W;
        Q = Q - X(:, i) * W;
    end
end