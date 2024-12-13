function [Sperp, np, S, B, Mind, sub_xkin] = choose_subspaces(X, delta, xkin, Pars, max_sketchsize)

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
            if np == max_sketchsize
                break  % Breaks out of for loop
            end
        end
    end
end

% initial guess of subspaces
S = Q(:, 1 : np);  
Sperp = Q(:, np + 1:n)';  % Will be empty if np=n

% special case
if isempty(S)
    Mind = xkin;
    sub_xkin = 1;
    B = [];
    return 
end

%% Now we have to check - do we have the points in X to guarantee good geometry
%% in this subspace? 

% Using the points projected onto the subspace, B, can I get a valid
% geometry? 
valid = 0;
vf = 1;
while ~valid
    [B, Mind, sub_xkin] = find_points_in_span(X, S, Sperp, xkin, delta, Pars, nf);
    [np, ~, ~, valid, ~, S_sub, Sperp_sub] = yet_another_formquad(B, [], delta, 0, Pars, sub_xkin, vf);
    if ~valid
        if isempty(S_sub)
            S = [];
            Sperp = eye(n);
            np = 0; 
            B = [];
            sub_xkin = 1; 
            Mind = xkin;
            valid = true;
            break
        else
            S = S * S_sub; 
            dim_sub = size(S, 2);
            % determine Sperp_sub
            [Q, ~] = qr(S);
            Sperp = Q(:, (dim_sub + 1):end)';
        end
    end
end


end