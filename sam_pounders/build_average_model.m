function [Cres, Gres, Hres] = build_average_model(models, X_inc, subset)

    % the subset argument exists for when you only want to compute and
    % store the data for just a subset of the models, which is a necessary
    % operation for computing the ameliorated model - you save the model
    % data about to be updated, update those models, and then get the full
    % average model. 

    if nargin == 2
        subset = 1:length(models);
    else
        subset = subset';
    end

    % fetch data from models
    m = length(subset);
    n = models(1).n;
    Cres = zeros(1, m); Gres = zeros(n, m); Hres = zeros(n, n, m);
    for j = subset
        Cj = models(j).Cres;
        Gj = models(j).Gres;
        Hj = models(j).Hres;
        center_point = models(j).center_point;

        % if this next step seems weird to you, just remember:
        % each model m_i has its own center_point c_i and so can be expressed
        % as
        % m_i(c_i + s) = C + G'*s + 0.5*s'*H*s
        % since the average model is centered at X_inc,
        % we want a quadratic model (using a variable d) to satisfy
        % m_i(X_inc + d) = m_i(c_i + s) after a change of variables.
        % That change of variables is thus:
        % s = X_inc + d - c_i = X_inc - c_i + d
        % and that's what we're applying here: 
        displaced = X_inc - center_point;
        Cres(j) = Cj + displaced * Gj + 0.5 * displaced * Hj * displaced';
        Gres(:, j) = Gj + (displaced * Hj)';
        Hres(:, :, j) = Hj;
    end
end