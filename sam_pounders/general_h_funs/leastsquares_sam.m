function [c, G, H] = leastsquares_sam(Cres, Gres, Hres, X, centers, xkin)
[n, ~, m] = size(Hres);

D = X(centers, :) - repmat(X(xkin, :), m, 1);

H = zeros(n); G = zeros(n,1); c = 0.5*sum(Cres.^2);
for i = 1:m
    g = 2 * Cres(i) * Gres(:, i);
    h = Cres(i) * Hres(:, :, i);
    H = H + 2 * h + 2 * (Gres(:, i) * Gres(:, i)');
    G = G + g;
    G = G - h*D(i, :)';
    c = c + 0.5 * D(i, :) * h * D(i, :)';
    c = c - D(i, :) * g;
end