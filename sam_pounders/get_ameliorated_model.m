function [G, H] = get_ameliorated_model(X, xkin, old_centers, oldCres, Cres, oldGres, Gres, Hres, Hresdel, models_to_update, probs, combinemodels)

    [n,~,m] = size(Hres);

    D = repmat(X(xkin, :), m, 1) - X(old_centers, :);

    G = zeros(n, 1); H = zeros(n);

    for i = 1:m
        if strcmp(functions(combinemodels).function,'leastsquares')
            old_ci = oldCres(i) + D(i, :) * oldGres(:, i) ...
                + 0.5 * D(i, :) * Hres(:, :, i) * D(i, :)';
            old_gi = old_ci*(oldGres(:, i) + Hres(:, :, i) * D(i, :)');
            old_Hi = old_ci*Hres(:, :, i);
            new_ci = Cres(i);
            new_gi = new_ci*Gres(:, i);
            new_Hi = new_ci*Hres(:, :, i) + new_ci*Hresdel(:, :, i);
            G = G + old_gi;
            H = H + old_Hi;
        elseif strcmp(functions(combinemodels).function,'identity_combine')
            old_gi = oldGres(:, i) + Hres(:, :, i) * D(i, :)';
            old_Hi = Hres(:, :, i);
            new_gi = Gres(:, i);
            new_Hi = Hres(:, :, i) + Hresdel(:, :, i);
            G = G + old_gi;
            H = H + old_Hi;
        end
        if ismember(i, models_to_update)
            G = G + (1.0/probs(i))*(new_gi - old_gi);
            H = H + (1.0/probs(i))*(new_Hi - old_Hi);
        end
    end
    
    if strcmp(functions(combinemodels).function,'leastsquares')
        weights = zeros(1,m);
        weights(models_to_update) = (1.0./probs(models_to_update)).^2;
        H = H + oldGres * diag(1.0 - weights) * oldGres' + Gres * diag(weights) * Gres';
    end
    

end