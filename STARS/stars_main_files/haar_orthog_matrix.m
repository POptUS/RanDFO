function [Z] = haar_orthog_matrix(dimension)
    A = randn(dimension);
    [Z, R] = qr(A);
    Z = Z * diag(sign(diag(R)));
end
