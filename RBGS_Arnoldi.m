function [Q, H] = RBGS_Arnoldi(A, v, s, p, Theta, basisFunc)
    %Inputs:
    %A          matrix of size (n, n)
    %v          start vector of size (n, 1)
    %s          step size
    %p          number of iterations, dimension of the final
    %           Krylov subspace is m = s * p + 1
    %Theta      random matrix of size (d, n)
    %basisFunc  Krylov subspace basis generating function
    %
    %Outputs:
    %Q          Theta-orthonormal matrix of size (n, m)
    %H          upper Hessenberg matrix of size (m, m-1)

    m = s * p + 1;
    n = size(A, 1);
    d = size(Theta, 1);
    V = zeros(n, m); %store basis vectors
    P = zeros(d, m); %store sketched basis vectors
    S = zeros(d, m); %store orthonormal sketched basis vectors
    R = zeros(m+1, m+1);
    Q = zeros (n, m);
    H = zeros (m, m-1);
    B = zeros (s+1, s, p);

    V(:, 1) = v;
    P(:, 1) = Theta * v;
    R(1, 1) = norm(P(:, 1));
    Q(:, 1) = V(:, 1) ./ R(1, 1);
    S(:, 1) = Theta * Q(:, 1);

    for j = 1:p
        i = s * (j-1) + 1;
        b = (i + 1) : (i + s);
        [V(:, b), B(:, :, j)] = basisFunc(A, Q(:, i), s);
        P(:, b) = Theta * V(:, b);
        R(1:i, b) = S(:, 1:i) \ P(:, b);
        Q(:, b) = V(:, b) - Q(:, 1:i) * R(1:i, b);
        S(:, b) = P(:, b) - S(:,1:i) * R(1:i, b);
        [~, R(b, b)] = qr(S(:, b), 0);
        Q(:, b) = Q(:, b) / R(b, b);
        S(:, b) = Theta * Q(:, b);
        b_0 = i : (i + s - 1);
        %update Hessenberg, explicit sketch
        M = Theta * A * Q(:, b_0);
        H(1:(i+s), b_0) = S(:, 1:(i+s)) \ M;
    end
end