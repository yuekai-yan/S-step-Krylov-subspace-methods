function [x, beta] = BGS_GMRES(A, s, p, basisFunc, b, ctol)
    %Inputs:
    %A          matrix of size (n, n)
    %v          start vector of size (n, 1)
    %s          step size
    %p          number of iterations, dimension of the final
    %           Krylov subspace is m = s * p + 1
    %basisFunc  Krylov subspace basis generating function
    %b          RHS of the A * x = b
    %ctol       convergence tolerence
    %
    %Outputs:
    %x          approximate solution of size (n, 1)
    %beta       store the norm of the residual

    m = s * p + 1;
    n = size(A, 1);
    %d = size(Theta, 1);
    V = zeros(n, m); %store basis vectors
    %P = zeros(d, m); %store sketched basis vectors
    %S = zeros(d, m); %store orthonormal sketched basis vectors
    R = zeros(m+1, m+1);
    Q = zeros (n, m);
    H = zeros (m, m-1); %store upper Hessenberg matrix
    B = zeros (s+1, s, p);

    x0 = zeros(n, 1);
    r0 = b - A * x0;
    %beta0 = norm(Theta * r0);
    beta0 = norm(r0);
    beta = beta0;
    e1 = zeros(m, 1);
    e1(1) = beta0;

    v = r0;
    V(:, 1) = v;
    %P(:, 1) = Theta * v;
    %R(1, 1) = norm(P(:, 1));
    R(1, 1) = norm(V(:, 1));
    Q(:, 1) = V(:, 1) ./ R(1, 1);
    %S(:, 1) = Theta * Q(:, 1);

    for j = 1:p
        i = s * (j-1) + 1;
        cols = (i + 1) : (i + s);
        [V(:, cols), B(:, :, j)] = basisFunc(A, Q(:, i), s);
        %P(:, cols) = Theta * V(:, cols);
        %R(1:i, cols) = S(:, 1:i) \ P(:, cols);
        R(1:i, cols) = Q(:, 1:i)' * V(:, cols);
        Q(:, cols) = V(:, cols) - Q(:, 1:i) * R(1:i, cols);
        %S(:, cols) = P(:, cols) - S(:,1:i) * R(1:i, cols);
        [Q(:, cols), R(cols, cols)] = qr(Q(:, cols), 0);
        %Q(:, cols) = Q(:, cols) / R(cols, cols); 
        %S(:, cols) = Theta * Q(:, cols);
        b_0 = i : (i + s - 1);

        %update Hessenberg, explicit sketch
        %M = Theta * A * Q(:, b_0);
        M = A * Q(:, b_0);
        H(1:(i+s), b_0) = Q(:, 1:(i+s))' * M;

        %solve the least-square problem
        y = H(1:(i+s), 1:(j*s)) \ e1(1:(i+s), 1);

        %compute the residual
        beta = [beta; ...
                norm(e1(1:(i+s), 1) - H(1:(i+s), 1:(j*s)) * y)];
        k = j * s;
        %test for convergence of residual
        if beta(end) < ctol * beta0
            fprintf('converged after %d steps\n', k);
            break
        end

    end

    %add correction from current Krylov subspace
    x = x0 + Q(:, 1:k) * y;
end