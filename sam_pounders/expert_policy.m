function probs = lipschitz_estimate_policy(models, batch_size, sample_type, data)

    X_inc = data.X_inc;
    delta = data.delta;
    already_updated = data.already_updated;

    m = length(models);
    error_estimates = zeros(1, m); 
    for j = 1:m
        if strcmp(sample_type, 'model')
            if ~ismember(j, already_updated)
                % we assign 0 probs to already_updated. 
                error_estimates(j) = get_error_estimate(models(j), x, delta);
            end
        elseif strcmp(sample_type, 'rho')
            error_estimate_inc = get_error_estimate(models(j), X_inc);
            error_estimate_trial = get_error_estimate(models(j), x); 
            error_estimates(j) = max(error_estimate_inc, error_estimate_trial);
        else
            error('sample_type input to expert_policy must be a string model or rho ')
        end
    end

    % optimal probabilities based on error estimates
    [probs, var] = compute_probs(batch_size, error_estimates);

end
