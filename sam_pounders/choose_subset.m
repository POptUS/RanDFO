function [probs, to_update, V] = choose_subset(models, batch_size, x, delta, already_updated)

    m = length(models);
    error_estimates = zeros(1, m); 
    for j = 1:m
        if nargin == 5
            error_estimates(j) = get_error_estimate(models(j), x, delta);
        else
            error_estimates(j) = get_error_estimate(models(j), x);
        end
    end

    % censor pointless-to-update models
    if nargin == 5
        error_estimates(already_updated) = 0;
    end

    % optimal probabilities based on error estimates
    [probs, V] = compute_probs(batch_size, error_estimates);

    % now correct the probabilities for Poisson sampling
    new_probs = get_cps_probs(probs', batch_size);
    if abs(sum(new_probs) - batch_size) > 1.0
        % THIS IS A SAFEGUARD IF SOMETHING WENT NUMERICALLY AWRY. IT
        % SHOULD TRIGGER RARELY.  
        % 1.0 IS ARBITRARY, BUT INDICATES SOMETHING HAS GONE QUITE
        % WRONG. 
        warning('Something went far too wrong in the CPS computation. Doing a greedy thing on this iteration.')
        [~, to_update] = maxk(new_probs, batch_size); 
        rejected = false;
    else
        rejected = true;
    end
    while rejected
        % POISSON SAMPLE
        to_update = [];
        for k = 1:m
            % flip coin
            if rand() < new_probs(k)
                to_update = cat(1,to_update,k);
            end
        end
        if length(to_update) == batch_size
            rejected = false;
        end
    end
    to_update = to_update(:);

end