function probs = surrogate_informed_policy(models, batch_size, sample_type, data, surrogate_array, Low, Upp, spsolver)

    X_inc = data.X_inc;
    if strcmp(sample_type, 'model')
        delta = data.delta;
        already_updated = data.already_updated;
    elseif strcmp(sample_type, 'rho')
        Xsp = data.Xsp;
    end

    m = length(models);
    reward_estimates = zeros(1, m); 
    for j = 1:m
        pre_Cres = models(j).Cres;
        pre_Gres = models(j).Gres;
        pre_Hres = models(j).Hres;
        center_point = models(j).center_point;
        if strcmp(sample_type, 'model')
            if ~ismember(j, already_updated)
                % we assign 0 probs to already_updated.
                displaced = X_inc - center_point;
                [~, post_Cres, post_Gres, post_Hres] = fake_update_model(models(j), @(x)predict(surrogate_array{j}, x), X_inc, delta);
                diff_model_Cres = post_Cres - pre_Cres - displaced * pre_Gres - 0.5 * displaced * pre_Hres * displaced';
                diff_model_Gres = post_Gres - pre_Gres - pre_Hres * displaced';
                diff_model_Hres = post_Hres - pre_Hres;
                reward_estimates(j) = max_abs_diff(diff_model_Cres, diff_model_Gres, diff_model_Hres, delta, Low, Upp, X_inc, spsolver, length(X_inc));
            end
        elseif strcmp(sample_type, 'rho')
            displaced = X_inc - center_point;
            reward_estimate_inc = abs(predict(surrogate_array{j}, X_inc) - pre_Cres - displaced * pre_Gres - 0.5 * displaced * pre_Hres * displaced');
            displaced = Xsp - center_point;
            reward_estimate_trial = abs(predict(surrogate_array{j}, Xsp) - pre_Cres - displaced * pre_Gres - 0.5 * displaced * pre_Hres * displaced');
            reward_estimates(j) = max(reward_estimate_inc, reward_estimate_trial);
        else
            error('sample_type input to surrogate_informed_policy must be a string model or rho ')
        end
    end

    % optimal probabilities based on error estimates
    [probs, var] = compute_probs(batch_size, reward_estimates);

end