function [X,F,have_eval,xkin] = condense(X,F,have_eval,nf_start,nf_end,xkin)

    to_remove = [];
    for j = nf_start:nf_end
        x = X(j, :);
        for i = 1:nf_end
            if all(isalmost(x, X(i, :), 1e-12))
                have_eval(i, :) = have_eval(j, :) + have_eval(i, :);
                F(i, F(j, :) ~= 0) = F(j, F(j, :) ~= 0);
                to_remove = cat(1, to_remove, j);
            end
        end
    end

    X(to_remove, :) = [];
    F(to_remove, :) = [];
    have_eval(to_remove, :) = [];
    
    xkin = xkin - sum(to_remove < xkin);

end