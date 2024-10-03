function H = vecToSymmetric(h, n)
    % vecToSymmetric reconstructs the symmetric matrix H from its vectorized upper triangular part h
    % Input:
    %   h - a vector containing the unique elements of the upper triangular part of H
    %   Therefore the dimension of h must be n(n+1)/2
    %   n - the dimension of the square matrix H
    % Output:
    %   H - the reconstructed n-by-n symmetric matrix

    H = zeros(n, n);
    idx = 1;

    % Fill the upper triangular part (including diagonal) of H
    for i = 1:n
        for j = i:n
            H(i, j) = h(idx);
            if i ~= j
                H(j, i) = h(idx);  % Ensure symmetry
            end
            idx = idx + 1;
        end
    end
end
