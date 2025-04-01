function [FX_inc, FX_sp] = evaluate_ameliorated_model(models, probs, to_update, X_inc, Xsp, old_F_inc, old_F_sp)

    m = length(models);
    FX_inc = zeros(1, m);
    FX_sp = zeros(1, m);

    for j = 1:m
        FX_inc(j) = model_value_at_point(models(j), X_inc);
        FX_sp(j) = model_value_at_point(models(j), Xsp);       
    end

    ctr = 1;
    for j = to_update'
        FX_inc(j) = (1.0/probs(j)) * FX_inc(j);
        FX_inc(j) = FX_inc(j) + (1.0 - 1.0/probs(j)) * old_F_inc(ctr);
        FX_sp(j) = (1.0/probs(j)) * FX_sp(j);
        FX_sp(j) = FX_sp(j) + (1.0 - 1.0/probs(j)) * old_F_sp(ctr);
        ctr = ctr + 1;
    end
end