function ind = check_if_evaluated(cand, X, nf_start, nf)

    ind = 0; 
    for j = nf_start:nf
        x = X(j, :);
        if all(isalmost(cand,x,eps))
            ind = j;
            break;
        end
    end

end