function [models, FX_inc, FX_sp, new_evals_ctr] = evaluate_ameliorated_model(models, probs, to_update, X_inc, Xsp, model_prediction_inc, model_prediction_Xsp, iter, delta, nf, nf_max)

    is_at_center = center_checker(models, X_inc);

    FX_inc = model_prediction_inc;
    FX_sp = model_prediction_Xsp;

    new_evals_ctr = 0;

    for j = to_update'
        % X_inc evaluation
        if ~is_at_center(j)
            % gotta update this model
            models(j) = update_center_point(models(j), X_inc);
            new_evals_ctr = new_evals_ctr + 1;
            models(j).trial_iters(iter) = 1;
            % update the model with new center point
            [models(j), new_evals] = update_model(models(j), delta, nf_max, nf + new_evals_ctr);
            new_evals_ctr = new_evals_ctr + new_evals;
        end
        Fj_inc = models(j).Cres;

        % Xsp evaluation
        models(j) = tentative_update_center_point(models(j), Xsp);
        Fj_sp = models(j).F(models(j).center_idx); 
        
        new_evals_ctr = new_evals_ctr + 1;
        if is_at_center(j)
            models(j).trial_iters(iter) = 1;
        else
            models(j).trial_iters(iter) = 2;
        end

        % populate FX_inc, FX_sp
        FX_inc(j) = (1.0 / probs(j)) * (Fj_inc - FX_inc(j)) + FX_inc(j);
        FX_sp(j) = (1.0 / probs(j)) * (Fj_sp - FX_sp(j)) + FX_sp(j);
    end
end