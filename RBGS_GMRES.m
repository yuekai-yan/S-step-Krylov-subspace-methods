function [x, beta, orthErr] = RBGS_GMRES(A, s, p, Theta, basisFunc, ...
                              AOB, WB, b, ctol)
    %Inputs:
    %A          matrix of size (n, n) or function handle
    %v          start vector of size (n, 1)
    %s          step size
    %p          number of iterations, dimension of the final
    %           Krylov subspace is m = s * p + 1
    %Theta      random matrix of size (d, n)
    %basisFunc  Krylov subspace basis generating function
    %AOB        orthogonalization methods against other blocks {
    %           RGS, rCGS, rMGS}
    %WB         orthogonalization methods within one block {
    %           rWhitening (QR), RGS (backslash), rCGS, rCGS2, rMGS}
    %b          RHS of the A * x = b
    %ctol       convergence tolerence
    %
    %Outputs:
    %x          approximate solution of size (n, 1)
    %beta       store the norm of the residual

    % check whether A is a matrix or function handle
    if isa(A, 'function_handle')
        Amul = @(x) A(x);
    else
        Amul = @(x) A * x;
    end

    m = s * p + 1;
    n = length(b);
    d = size(Theta, 1);
    V = zeros(n, m); % store basis vectors
    %P = zeros(d, m); % store sketched basis vectors
    S = zeros(d, m); % store orthonormal sketched basis vectors
    R = zeros(m+1, m+1);
    Q = zeros (n, m);
    H = zeros (m, m-1);
    B = zeros (s+1, s, p);
    orthErr = 0; % store the orthogonalization error of the 
                  % Theta-orthogonal matrix Q
                  
    x0 = zeros(n, 1);
    r0 = b - Amul(x0);
    beta0 = norm(Theta * r0);
    beta = norm(r0);
    e1 = zeros(m, 1);
    e1(1) = beta0;

    v = b;
    V(:, 1) = v;
    sketch0 = Theta * v;
    R(1, 1) = norm(sketch0);
    Q(:, 1) = V(:, 1) ./ R(1, 1);
    S(:, 1) = Theta * Q(:, 1);

    for j = 1:p
        i = s * (j-1) + 1;
        cols = (i + 1) : (i + s);
        k = j * s;
        [V(:, cols), B(:, :, j)] = basisFunc(Amul, Q(:, i), s);

        % Theta-orthogonalized V(:, cols) with Q(:, 1:i)
        AOB_fun = str2func(sprintf('AOB.%s', AOB));
        [Q(:, cols), R(1:i, cols)] = AOB_fun(Q(:, 1:i), ...
            V(:, cols), Theta);

        % Theta orthonormal Q(:, cols)
        WB_fun  = str2func(sprintf('WB.%s', WB));
        [Q(:, cols), R(cols, cols)] = WB_fun(Q(:, cols), Theta);

        Q1 = Theta * Q(:, 1:k);
        orthErr = [orthErr; norm(Q1' * Q1 - eye(k), 'fro')];
        S(:, cols) = Theta * Q(:, cols);
        b_0 = i : (i + s - 1);  % s*(j-1)+1 : s*j
        
        % update Hessenberg, explicit version
        M = Theta * Amul(Q(:, b_0));
        H(1:(i+s), b_0) = S(:, 1:(i+s)) \ M;

        % solve the least-square problem
        y = H(1:(i+s), 1:(j*s)) \ e1(1:(i+s), 1);  % k-by-1

        % compute the residual
        beta = [beta; ...
                norm(Amul(Q(:,1:length(y))*y)-b) / norm(b)];
        % test for convergence of residual
        if beta(end) < ctol
            fprintf('RBGS converged after %d steps\n', k);
            break
        end

    end

    % add correction from current Krylov subspace
    x = x0 + Q(:, 1:k) * y;
end