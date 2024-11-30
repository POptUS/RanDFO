function M = matYYT(Y)
    % matYYT computes the vectorized upper triangular part of yy^T, where y
    % is a data point in Y
    % Input:
    % Y - a m-by-n matrix of m data points of dimension n: Y = [y1; y2; ...; ym]
    % where yi is an n-dimensional vector, for i = 1, 2, ..., m
    % Output:
    % The ith row of M is a vector containing the unique elements of the upper triangular part of yi * yi^T
    % with crossed terms in yi * yi^T multiplied by 2

    [m, n] = size(Y);
    for k = 1:m
        Z(:, :, k) = Y(k, :)' * Y(k, :);
    end

    idx = 1;
    for i = 1:n
        for j = i:n
            if i == j
                M(:, idx) = Z(i, j, :);
            else
                M(:, idx) = 2 * Z(i, j, :);
            end
            idx = idx + 1;
        end
    end
end
