function subset = computeUCBsubset(mean_bandit, sketch_dictionary, regularizer, Mdir, coeff, threshold, max_sketchsize)

    % special case handling:
    coeff = max(eps, coeff);

    dim = size(Mdir, 1);
    % special case handling:
    if dim == 0
        subset = [];
        return
    end

    scores = zeros(1, dim);
    perm = randperm(dim, min(dim, max_sketchsize * 10));

    for k = perm %1:dim
        q = Mdir(k, :);
        icq = smw_product(sketch_dictionary, regularizer, q');
        scores(k) = abs(q * mean_bandit) + coeff * sqrt(q * icq);
    end

    %sorted_percents = cumsum(sort(scores, 'descend')) / sum(scores);
    %k = find(sorted_percents > threshold, 1);
    %[~, subset] = maxk(scores, min(k, max_sketchsize));
    
    [~, subset] = maxk(scores, max_sketchsize);

end
