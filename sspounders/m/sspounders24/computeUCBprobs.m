function pi = computeUCBprobs(mean_bandit, ic, Mdir, ss, coeff, barG)

    dim = size(Mdir, 1);

    scores = zeros(1, dim);

    for k = 1:dim
        q = Mdir(k, :);
        scores(k) = max(abs(q * (barG - mean_bandit) + coeff * sqrt(q * ic * q')), ...
            abs(q * (barG - mean_bandit) - coeff * sqrt(q * ic * q')));
    end

    scores = max(scores.^2, eps);

    [sorted_scores, sort_inds] = sort(scores); 
    cum_sorted_scores = cumsum(sorted_scores);

    % find largest c satisfying lhs <= rhs
    lhs = (1:p) + ss - p;
    rhs = cum_sorted_scores./sorted_scores;
    c = find((lhs<=rhs).*(lhs>0),1,'last');
    
    % now compute the probabilities
    pi = ones(p,1);
    if ~isempty(c)
        pi(sort_inds(1:c)) = (lhs(c)/cum_sorted_scores(c))*sorted_scores(1:c);
    end
    
    
    % to avoid divide by zero issues
    pi(pi <= eps) = eps;
end