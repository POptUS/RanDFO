function [Cres, Gres, Hres] = build_ameliorated_model(models, X_inc, subset, probs, Cres, Gres, Hres)

    % input Cres, Gres, and Hres are holding the average model data from
    % before any models were updated

    % by the time this function is called, the models have been updated, so
    % this is computing m_i(x; X_inc) 
    [new_Cres, new_Gres, new_Hres] = build_average_model(models, X_inc, subset);

    batch_size = length(subset);
    for j = 1:batch_size
        Cres(subset(j)) = (1.0 - 1.0/probs(subset(j))) * Cres(subset(j));
        Cres(subset(j)) = Cres(subset(j)) + (1.0/probs(subset(j))) * new_Cres(j);
        Gres(:, subset(j)) = (1.0 - 1.0/probs(subset(j))) * Gres(:, subset(j));
        Gres(:, subset(j)) = Gres(:, subset(j)) + (1.0/probs(subset(j))) * new_Gres(:, j);
        Hres(:, :, subset(j)) = (1.0 - 1.0/probs(subset(j))) * Hres(:, :, subset(j));
        Hres(:, :, subset(j)) = Hres(:, :, subset(j)) + (1.0/probs(subset(j))) * new_Hres(:, :, j);
    end
end