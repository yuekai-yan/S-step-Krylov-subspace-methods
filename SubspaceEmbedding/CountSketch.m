function X = CountSketch(d, n)
    rowIndex = randi(d, 1, n);
    signs = randi(2, 1, n) * 2 - 3;
    X = sparse(rowIndex, 1:n, signs, d, n);
end