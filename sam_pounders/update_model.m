function [flag, X, F, nf, ncf, Gres, Hresdel, have_eval, valid, freebies, ncf_vec] = ...
    update_model(fun, X, F, L, U, delta, xkin, nf, ncf, Cres, Gres, Hres, models_to_update, have_eval, centers, npmax, nfmax, Par, valid, nf_start, must_update, ncf_vec)

[n, ~, m] = size(Hres);

flag = []; 

%if ~must_update
    freebies = [];
%end

% 1a. Compute the interpolation set for the models to update.
    Hresdel = zeros(n, n, m);
    
    % Need to do this for each model to update. 
    for i = models_to_update'
        Xi = X(have_eval(1:nf,i), :);
        Fi = F(have_eval(1:nf,i), i);
        numx = size(Xi, 1);
        Resi = zeros(numx, 1);
        row = 1:nf; xkini = find(row(have_eval(1:nf,i)) == centers(i)); 

        for j = 1:numx
            D = Xi(j, :) - Xi(xkini, :);
            Resi(j) = Fi(j) - Cres(i) - .5 * D * Hres(:, :, i) * D';
        end
         
        [Mdir, np, validi, Gresi, Hresdeli, Mind] = ...
            formquad(Xi, Resi, delta, xkini, npmax, Par, 0);
        valid(i) = validi;
        
        if np >= n %validi
            Gres(:, i) = Gresi; Hresdel(:, :, i) = Hresdeli;
        end
        if must_update
            if np < n  % Must obtain and evaluate bounded geometry points
                [Mdir, np] = bmpts(X(xkin, :), Mdir(1:n - np, :), L, U, delta, Par(3));
                for j = 1:min(n - np, nfmax - nf)
                    cand = min(U, max(L, X(xkin, :) + Mdir(j, :)));
                    ind = check_if_evaluated(cand, X, nf_start, nf);
                    if ind == 0
                        nf = nf + 1;
                        X(nf, :) = min(U, max(L, X(xkin, :) + Mdir(j, :))); % Temp safeguard
                        ind = nf;
                    end
                    Xi = cat(1, Xi, X(ind, :));
                    if ~have_eval(ind, i) 
                        F(ind, i) = fun(X(ind, :), i);
                        ncf = ncf + 1;
                        ncf_vec(nf) = ncf;
                        have_eval(ind, i) = true; 
                    end
                    if any(isnan(F(ind, :)))
                        [X, F, flag] = prepare_outputs_before_return(X, F, nf, -3);
                        return
                    end
                    D = Mdir(j, :);
                    Res_nf_i = F(ind, i) - Cres(i) - .5 * D * Hres(:, :, i) * D';
                    Resi = cat(1, Resi, Res_nf_i);
                end
                if ncf >= nfmax
                    break
                end
                [~, np, validi, Gresi, Hresdeli, Mind] = ...
                    formquad(Xi, Resi, delta, xkini, npmax, Par, 0);
                Gres(:, i) = Gresi; Hresdel(:, :, i) = Hresdeli;
                
                valid(i) = validi; 
                if np < n
                    [X, F, flag] = prepare_outputs_before_return(X, F, nf, -5);
                    return
                end
            end
        else % not must_update
            if np >= n
                freebies = cat(1,freebies,i);
            end
        end
    end
end