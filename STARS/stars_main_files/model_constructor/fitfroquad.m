function [G, H] = fitfroquad(poised_set, f_poised)
% These assumed that you are inputing at least n+1 poised points.
warning('off', 'all');
n = size(poised_set, 2);

M = [ones(n + 1, 1) poised_set(1:n + 1, :)]';

if size(poised_set, 1) == n + 1
    H = zeros(n, n);
    B = M'\f_poised(:)
    G = B(2:n + 1);

elseif size(poised_set, 1) > n+1
    N = zeros(.5 * n * (n + 1), n + 1);
    for np = 1:n + 1
        N(:, np) = phi2eval(poised_set(np, :))';
    end

    [Q, R] = qr(M');

    for i = n + 2:size(poised_set, 1)

        Ny = [N phi2eval(poised_set(i, :))'];
        [Qy, Ry] = qrinsert(Q, R, np + 1, [1 poised_set(i, :)], 'row'); % Update QR
        Ly = Ny * Qy(:, n + 2:np + 1);

        np = np + 1;
        N = Ny;
        Q = Qy;
        R = Ry;
        L = Ly;

        Z = Q(:, n + 2:np);
        M = [M [1; poised_set(i, :)']]; % Note that M is growing
    end

    % For L=N*Z, solve L'*L*Omega = Z'*f_poised:
    Omega = L' \ (Z' * f_poised(:));
    Omega = L \ Omega;
    Beta = L * Omega;
    Alpha = M' \ (f_poised(:) - N' * Beta);

    G(:) = Alpha(2:n + 1);
    num = 0;
    for i = 1:n
        num = num + 1;
        H(i, i) = Beta(num);
        for j = i + 1:n
            num = num + 1;
            H(i, j) = Beta(num) / sqrt(2);
            H(j, i) = H(i, j);
        end
    end
end

    function Phi = phi2eval(X)
        [m, n] = size(X);
        Phi = zeros(m, .5 * n * (n + 1));
        j = 0;
        for k = 1:n
            j = j + 1;
            Phi(:, j) = .5 * X(:, k).^2;
            for kk = k + 1:n
                j = j + 1;
                Phi(:, j) = X(:, k) .* X(:, kk) / sqrt(2);
            end
        end
    end

end
