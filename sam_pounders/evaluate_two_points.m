function [models, FX_inc, FXsp, average_values_inc, average_values_sp, new_evals_ctr] = evaluate_two_points(models, X_inc, Xsp, delta, to_update, probs, iter)

    m = length(models);

    is_at_center = center_checker(models, X_inc);

    average_values_inc_old = zeros(1, m); average_values_sp_old = zeros(1, m);
    average_values_inc = zeros(1, m); average_values_sp = zeros(1, m);
    FX_inc = zeros(1, m); FXsp = zeros(1, m);
    new_evals_ctr = 0;

    for j = 1:m
        % first populate current values
        average_values_inc_old(j) = model_value_at_point(models(j), X_inc);
        average_values_sp_old(j) = model_value_at_point(models(j), Xsp);
    end
     

    for j = to_update'
        % evaluation at the center
        if is_at_center(j)
            FX_inc(j) = models(j).F(models(j).center_idx);
        else
            models(j) = tentative_update_center_point(models(j), X_inc);
            FX_inc(j) = models(j).F(models(j).center_idx);
            new_evals_ctr = new_evals_ctr + 1;
            % A tricky thing: If we can update the model center without
            % incurring additional function evaluations, we should just do
            % it! 
            [acceptable, ~] = probe_interpolation(models(j), delta);
            if acceptable
                [models(j), new_evals] = update_model(models(j), delta, Inf, 0); % last two args should do nothing. 
                new_evals_ctr = new_evals_ctr + new_evals; % new_evals_should be 0, but let's just be safe! 
            else
                models(j) = make_untentative(models(j));
                % snap the model back to its previous center
                [models(j), new_evals] = update_model(models(j), delta, Inf, 0);
                new_evals_ctr = new_evals_ctr + new_evals; % new_evals should actually be 0 here, but just in case!
            end
        end
        % evaluation at the trial
        models(j) = tentative_update_center_point(models(j), Xsp);
        FXsp(j) = models(j).F(models(j).center_idx);
        models(j).trial_iters(iter) = 1;
        new_evals_ctr = new_evals_ctr + 1;
    end

    not_to_update = setdiff(1:m, to_update);

    for j = not_to_update
        FX_inc(j) = model_value_at_point(models(j), X_inc);
        FXsp(j) = model_value_at_point(models(j), Xsp);
    end
    FX_inc(not_to_update) = average_values_inc(not_to_update);
    FXsp(not_to_update) = average_values_sp(not_to_update); 

    % return average values
    average_values_inc(to_update) = FX_inc(to_update);
    average_values_inc(not_to_update) = average_values_inc_old(not_to_update);
    average_values_sp(to_update) = FXsp(to_update); 
    average_values_sp(not_to_update) = average_values_sp_old(not_to_update);

    % if we are ameliorating the estimates:
    FX_inc(to_update) = (1.0./probs(to_update)').*FX_inc(to_update) + (1.0 - 1.0./probs(to_update)').*average_values_inc_old(to_update);
    FXsp(to_update) = (1.0./probs(to_update)').*FXsp(to_update) + (1.0 - 1.0./probs(to_update)').*average_values_sp_old(to_update);

end