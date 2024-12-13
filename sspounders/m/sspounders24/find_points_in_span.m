%routine to find all points in X \intersect {x_k + S's: s\in R^(sub_dim)}

% returns 
% Mind: indices of X, Mind satisfying this property and
% 

function [B, Mind, sub_xk_in] = find_points_in_span(X, S, Sperp, xkin, delta, Pars, nf)

    if isempty(Sperp)
        Mind = 1:nf;
        B = (X(1: nf, :) - repmat(X(xkin, :), nf, 1)) * S;
        sub_xk_in = xkin;
        return
    end

    % special case
    if isempty(S)
        Mind = xkin;
        B = zeros(1, size(Sperp, 1));
        sub_xk_in = 1;
        return
    end
    
    dim = size(S, 2);
    B = zeros(1, dim);
    np = 0;
    Mind = [];
    
    for j = 1:nf 
        d = X(j, :) - X(xkin, :);
        if (norm(d) < delta * Pars(2) && norm(Sperp * d') < delta * sqrt(eps))
            np = np + 1;
            B(np, :) = (d * S)';
            Mind = [Mind, j];                     
        end
        if j == xkin
            sub_xk_in = length(Mind);
        end
    end
