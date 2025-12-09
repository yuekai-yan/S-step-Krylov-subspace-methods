function Afun = getStartMatrix(n, alpha, kappa)
%   [Afun, info] = getStartMatrix(n, kappa, alpha)
%   - n     : problem size
%   - kappa : 2-norm condition number
%   - alpha : weight of the upper diagonal (default = 0.1)
%
%   Returns:
%     Afun  : function handle for y = A*x with A = Q^H * B * Q, Q = FFT matrix
%     info  : struct with fields B (upper bidiagonal), kind, kappa_theory
%
%   Notes:
%   * Q is unitary (FFT), so cond2(A) = cond2(B).
%   * If alpha ≠ 0, the operator is generally non-symmetric (non-normal).

    if nargin < 2
        alpha = 0;
    end  
    if nargin < 3
        kappa = -3;  % default superdiagonal weight
    end

    % Construct eigenvalues (positive, logarithmically spaced)
    % Range: [1, kappa], length = n
    %lambda = linspace(1, kappa, n).';
    lambda = logspace(0,kappa,n).';

    % Build upper bidiagonal matrix in Fourier domain
    % B = diag(lambda) + alpha * diag(ones(n-1,1), 1);
    B = diag(lambda) + alpha * diag(randn(n-1,1), 1);
    % B = diag(lambda);

    % Define matvec: A*x = Q^H * B * Q * x
    Afun = @(x) ifft(B * fft(x));

    % Additional information
    %info = struct();
    %info.kind = 'Upper bidiagonal via unitary similarity: A = Q^H * B * Q';
    %info.B = B;
    %info.kappa_theory = cond(B);   % numerical condition number
end