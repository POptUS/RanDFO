function Q = compute_Q_factor(X, delta, xkin, Pars)

% Internal parameters:
[nf,n] = size(X);

% Precompute the scaled displacements
D = zeros(nf, n); % Stores scaled displacements
Nd = zeros(nf, 1); % Stores norm of scaled displacements
for i = 1:nf
    D(i, :) = (X(i, :) - X(xkin, :)) / delta;
    Nd(i) = norm(D(i, :));
end

% Try to find n+1 sufficiently affinely independent points:
Q = eye(n);
R = []; % Initialize the QR factorization of interest
Mind = xkin; % Indices of model interpolation points
np = 0;  % Counter for number of interpolation points

for i = nf:-1:1
    if Nd(i) <= Pars(1)
        proj = norm(D(i, :) * Q(:, np + 1:n), 2); % Project D onto null
        if proj >= Pars(3)  % add this index to Mind
            np = np + 1;
            Mind(np + 1, 1) = i;
            % MATLAB:
            % [Q,R] = qrinsert(Q,R,np,D(i,:)'); % Update QR
            % This bit is just for Octave
            if size(R, 1) == 0
                [Q, R] = qr(D(i, :)');
            else
                [Q, R] = qrinsert(Q, R, np, D(i, :)'); % Update QR
            end
            if np == n
                break  % Breaks out of for loop
            end
        end
    end
end

end