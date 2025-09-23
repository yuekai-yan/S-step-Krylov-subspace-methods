function X = Gaussian(d, n, mu, sigma)
    %X is a d-by-n random matrix with independent entries
    %distributed as N(mu,sigma^2).

    if nargin < 3
        mu = 0;
    end
    
    if nargin < 4
        sigma = 1;
    end

    X = mu + sigma .* randn(d, n);
end
