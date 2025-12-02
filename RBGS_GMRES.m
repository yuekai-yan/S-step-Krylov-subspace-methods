function [x, relRes, orthLoss, cond_num] = RBGS_GMRES(A, s, p, Theta, basisFunc, ...
                              AOB, WB, b, ctol)
    %Inputs:
    %A          matrix of size (n, n) or function handle
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
    %relRes       store the norm of the relative residual
    %orthLoss      store the loss of orthogonality
    %cond_num      store the condition number of Q

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
    Q = zeros(n, m);
    H = zeros(m, m-1);
    B = zeros(s+1, s);
    %M = [];  % store the Theta * A * Q
    orthLoss = 0; % store the orthogonalization error of the 
                  % Theta-orthogonal matrix Q
                  
    x0 = zeros(n, 1);
    x = x0;
    r0 = b - Amul(x0);
    beta0 = norm(Theta * r0);
    relRes = norm(r0) / norm(b);
    e1 = zeros(m, 1);
    e1(1) = beta0;

    v = b;
    V(:, 1) = v;
    sketch0 = Theta * v;
    %P(:, 1) = sketch0;
    R(1, 1) = norm(sketch0);
    Q(:, 1) = V(:, 1) ./ R(1, 1);
    cond_num = 1;
    S(:, 1) = Theta * Q(:, 1);

    for j = 1:p
        i = s * (j-1) + 1;
        cols = (i + 1) : (i + s);
        k = j * s;
        [V(:, cols), B] = basisFunc(Amul, Q(:, i), s);
        %P(:, cols) = Theta * V(:, cols);

        % Theta-orthogonalized V(:, cols) with Q(:, 1:i)
        AOB_fun = str2func(sprintf('AOB.%s', AOB));
        [Q(:, cols), R(1:i, cols)] = AOB_fun(Q(:, 1:i), ...
            S(:, 1:i), V(:, cols), Theta);

        % Theta orthonormal Q(:, cols)
        WB_fun  = str2func(sprintf('WB.%s', WB));
        [Q(:, cols), R(cols, cols)] = WB_fun(Q(:, cols), Theta);
        
        S(:, cols) = Theta * Q(:, cols);
        cond_num = [cond_num; cond(Q(:, 1:(k+1)))];
        orthLoss = [orthLoss; norm(S(:, 1:(k+1))' * S(:, 1:(k+1)) - eye(k+1), 'fro')];
        b_0 = i : k;  % s*(j-1)+1 : s*j
        
        % update Hessenberg, explicit version
        M = Theta * Amul(Q(:, b_0));
        H(1:(i+s), b_0) = S(:, 1:(i+s)) \ M;
        %{
        % implicit version
        b_hat = (i+1) : (i+s-1);
        K0 = [S(:, i), P(:, cols)] * B;
        M = [M, K0(:, 1)];
        K = K0(:, 2:end) - M(:, 1:i) * R(1:i, b_hat);
        K1 = K / R(b_hat, b_hat);
        M = [M, K1];
        H(1:(i+s), b_0) = S(:, 1:(i+s)) \ M(:, b_0); 
        %}

        % solve the least-square problem
        y = H(1:(i+s), 1:(j*s)) \ e1(1:(i+s), 1);  % k-by-1

        % compute the residual
        relRes = [relRes; ...
                norm(Amul(Q(:,1:length(y))*y)-b) / norm(b)];

        % add correction from current Krylov subspace
        x = [x, x0 + Q(:, 1:k) * y];

        % test for convergence of residual
        if relRes(end) < ctol
            fprintf('RBGS converged after %d steps\n', k);
            break
        end

    end

    % add correction from current Krylov subspace
    %x = x0 + Q(:, 1:k) * y;
end