function probs = lipschitz_estimate_policy(models, batch_size, sample_type, data)

    X_inc = data.X_inc;
    hfun = data.hfun;
    if strcmp(sample_type, 'model')
        delta = data.delta;
        already_updated = data.already_updated;
    elseif strcmp(sample_type, 'rho')
        Xsp = data.Xsp;
    end

    m = length(models);
    error_estimates = zeros(1, m); 
    for j = 1:m
        if strcmp(sample_type, 'model')
            if ~ismember(j, already_updated)
                % we assign 0 probs to already_updated. 
                error_estimates(j) = get_error_estimate(models(j), X_inc, delta);
            end
        elseif strcmp(sample_type, 'rho')
            error_estimate_inc = get_error_estimate(models(j), X_inc);
            error_estimate_trial = get_error_estimate(models(j), Xsp); 
            error_estimates(j) = max(error_estimate_inc, error_estimate_trial);
        else
            error('sample_type input to expert_policy must be a string model or rho ')
        end
        if strcmp(func2str(data.hfun), 'leastsquares')
            error_estimates(j) = error_estimates(j) * models(j).Cres;
        end
    end

    % optimal probabilities based on error estimates
    [probs, var] = compute_probs(batch_size, error_estimates);

end
