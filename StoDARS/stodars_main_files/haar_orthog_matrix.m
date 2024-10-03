function [Z] = haar_orthog_matrix(n)
    A = randn(n);
    [Z, R] = qr(A);
    Z = Z * diag(sign(diag(R)));
end
