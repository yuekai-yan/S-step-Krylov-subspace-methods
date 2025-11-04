function [x, beta, orthErr] = BGS_GMRES0(A, s, p, basisFunc, ...
                                        AOB, WB, b, ctol)
    %Inputs:
    %A          matrix of size (n, n) or function handle
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
    %beta       relative residual || A * x - b || / || b ||

    % check whether A is a matrix or function handle
    if isa(A, 'function_handle')
        Amul = @(x) A(x);
    else
        Amul = @(x) A * x;
    end

    m = s * p + 1;
    n = length(b);
    %d = size(Theta, 1);
    V = zeros(n, m);  % store basis vectors
    R = zeros(m+1, m+1);
    Q = zeros (n, m);
    H = zeros (m, m-1);  % store upper Hessenberg matrix
    B = zeros (s+1, s, p);

    x0 = zeros(n, 1);
    r0 = b - Amul(x0);
    beta0 = norm(r0);
    beta = 1;
    e1 = zeros(m, 1);
    e1(1) = beta0;

    orthErr = 0;  % store the orthogonalization error of the 
                  % Theta-orthogonal matrix Q

    v = r0;
    V(:, 1) = v;
    R(1, 1) = norm(V(:, 1));
    Q(:, 1) = V(:, 1) ./ R(1, 1);

    for j = 1:p
        i = s * (j-1) + 1;
        cols = (i + 1) : (i + s);
        k = j * s;
        [V(:, cols), B(:, :, j)] = basisFunc(Amul, Q(:, i), s);


        % orthogonalize V(:, cols) with Q(:, 1:i)
        AOB_fun = str2func(sprintf('AOB.%s', AOB));
        [Q(:, cols), R(1:i, cols)] = AOB_fun(Q(:, 1:i), ...
            V(:, cols));

        % orthonormalize Q(:, cols)
        WB_fun  = str2func(sprintf('WB.%s', WB));
        [Q(:, cols), R(cols, cols)] = WB_fun(Q(:, cols));

        

        %R(1:i, cols) = Q(:, 1:i)' * V(:, cols);
        %Q(:, cols) = V(:, cols) - Q(:, 1:i) * R(1:i, cols);

        % reorthogonalization between blocks
        %D = Q(:, 1:i)' * Q(:, cols);
        %Q(:, cols) = Q(:, cols) - Q(:, 1:i) * D;
        %R(1:i, cols) = R(1:i, cols) + D;

        % orthogonalization within a block
        %[Q(:, cols), R(cols, cols)] = qr(Q(:, cols), 0);

        % reorthogonalization within one block
        %[Q(:, cols), R1] = qr(Q(:, cols), 0);
        %R(cols, cols) = R1 * R(cols, cols);

        orthErr = [orthErr; norm(Q(:, 1:k)' * Q(:, 1:k) - eye(k), 'fro')];
        b_0 = i : (i + s - 1);

        % update Hessenberg, explicit sketch
        M = Amul(Q(:, b_0));
        H(1:(i+s), b_0) = Q(:, 1:(i+s))' * M;

        % solve the least-square problem
        y = H(1:(i+s), 1:(j*s)) \ e1(1:(i+s), 1);

        % compute the residual
        beta = [beta; ...
                norm(Amul(Q(:, 1:k)*y)-b) / norm(b)];
        % test for convergence of residual
        if beta(end) < ctol
            fprintf('BGS converged after %d steps\n', k);
            break
        end

    end

    % add correction from current Krylov subspace
    x = x0 + Q(:, 1:k) * y;
end