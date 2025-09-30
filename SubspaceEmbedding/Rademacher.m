function X = Rademacher(d, n)
    %X is a d-by-n matrix with independent entries
    %equal to +/- 1/sqrt(d) with probabilities 1/2
    R = randi([0, 1], d, n) * 2 - 1;
    X = R / sqrt(d);
end