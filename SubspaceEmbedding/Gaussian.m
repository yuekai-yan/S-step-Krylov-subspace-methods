function X = Gaussian(d, n)
    %X is a d-by-n random matrix with independent entries
    %distributed as N(0, 1/d).
    mu = 0;
    sigma = 1 / sqrt(d);
    %sigma = 1;
    X = mu + sigma .* randn(d, n);
end
