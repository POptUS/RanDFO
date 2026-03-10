function surrogate_array = rbf_surrogates(handle, n, m, lb, ub)

    % handle is for function you want to evaluate
    % n is the dimension of the input to handle
    % m is the dimension of the output of handle
    % lb is the row vector of lower bounds on the input
    % ub is the row vector of upper bounds on the input

    num_points = 500;

    X_unit = lhsdesign(num_points, n);

    X = lb + X_unit .* (ub - lb);

    % evaluate objective function
    F = zeros(num_points, m); 
    for i = 1:num_points
        F(i, :) = handle(X(i, :));
    end

    % build surrogates
    surrogate_array = cell(1, m);
    for j = 1:m 
        surrogate_array{j} = fitrgp(X, F(:, j), 'KernelFunction', 'squaredexponential');
    end
end