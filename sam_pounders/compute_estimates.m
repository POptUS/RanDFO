function [f0, fs] = compute_estimates(X, xkin, nf, Cres, Fxkin, Fy, centers, Gres, Hres, fc_to_update, fc_probs, fs_to_update, fs_probs, combinemodels)

    m = length(fc_probs);

    f0 = 0; fs = 0;

    Dc = repmat(X(xkin, :), m, 1) - X(centers, :);
    Ds = repmat(X(nf, :), m, 1) - X(centers, :);

    for j = 1:m
        old_cj = Cres(j) + Dc(j, :) * Gres(:, j) + 0.5 * Dc(j, :) * Hres(:, :, j) * Dc(j, :)';
        if strcmp(functions(combinemodels).function,'leastsquares')
            f0 = f0 + old_cj^2; 
            if ismember(j, fc_to_update)
                f0 = f0 + (1.0/fc_probs(j))*(Fxkin(j)^2 - old_cj^2);
                %f0 = f0 + (Fxkin(j)^2 - old_cj^2);
            end
        elseif strcmp(functions(combinemodels).function,'identity_combine')
            f0 = f0 + old_cj;
            if ismember(j, fc_to_update)
                f0 = f0 + (1.0/fc_probs(j))*(Fxkin(j) - old_cj);
                %f0 = f0 + (Fxkin(j) - old_cj);
            end
        end
        old_cj = Cres(j) + Ds(j, :) * Gres(:, j) + 0.5 * Ds(j, :) * Hres(:, :, j) * Ds(j, :)';
        if strcmp(functions(combinemodels).function,'leastsquares')
            fs = fs + old_cj^2;
            if ismember(j, fs_to_update)
                fs = fs + (1.0/fs_probs(j))*(Fy(j)^2 - old_cj^2);
                %fs = fs + (Fy(j)^2 - old_cj^2);
            end
        elseif strcmp(functions(combinemodels).function,'identity_combine')
            fs = fs + old_cj;
            if ismember(j, fs_to_update)
                fs = fs + (1.0/fs_probs(j))*(Fy(j) - old_cj);
                %fs = fs + (Fy(j) - old_cj);
            end
        end
    end

end