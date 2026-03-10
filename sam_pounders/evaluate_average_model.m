function [models, FX_inc, FXsp, new_evals_ctr] = evaluate_average_model(models, to_update, X_inc, Xsp, delta, nf_max, nf, iter)

    m = length(models);

    is_at_center = center_checker(models, X_inc);

    FX_inc = zeros(1, m); FXsp = zeros(1, m);
    new_evals_ctr = 0;

    for j = 1:m
        % evaluation at the center
        if is_at_center(j)
            FX_inc(j) = models(j).Cres;
        else
            if ismember(j, to_update)
                % need to evaluate! 
                % update the center
                models(j) = update_center_point(models(j), X_inc);
                new_evals_ctr = new_evals_ctr + 1;
                models(j).trial_iters(iter) = 1;
                % update the model with new center point
                [models(j), new_evals] = update_model(models(j), delta, nf_max, nf + new_evals_ctr);
                new_evals_ctr = new_evals_ctr + new_evals;
                FX_inc(j) = models(j).Cres;
            else
                % model fill-in
                FX_inc(j) = model_value_at_point(models(j), X_inc);
            end
        end
        % evaluation at the trial
        if ismember(j, to_update)
            % evaluate for real, make tentative center
            models(j) = tentative_update_center_point(models(j), Xsp);
            FXsp(j) = models(j).F(models(j).center_idx); 
            new_evals_ctr = new_evals_ctr + 1;
            if is_at_center(j)
                models(j).trial_iters(iter) = 1;
            else
                models(j).trial_iters(iter) = 2;
            end
        else
            % model fill-in
            FXsp(j) = model_value_at_point(models(j), Xsp);
        end
    end

end