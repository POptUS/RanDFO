function [new_probs, to_update] = poisson_sample_shortcut(probs, batch_size)

    % This is NOT technically accurate. This is a heuristic shortcut to
    % poisson_sample, because it is an extreme bottleneck in running this
    % code. 

    m = length(probs);
    new_probs = probs';

    % normalize
    new_probs = batch_size * new_probs / sum(new_probs); 

    if batch_size == m
        to_update = (1:m);
    else
        to_update = [];
        coinflips = zeros(1, m);
        for k = 1:m
            % flip coin
            coinflips(k) = rand();
            if coinflips(k) < new_probs(k)
                to_update = cat(1, to_update, k);
            end
        end
        len_to_update = length(to_update);
        if len_to_update > batch_size
            excess = len_to_update - batch_size; 
            [~, sorted_coinflips] = sort(coinflips);
            removed = 0; ctr = 1;
            while removed < excess
                if ismember(sorted_coinflips(ctr), to_update)
                    to_update(find(to_update == sorted_coinflips(ctr))) = [];
                    removed = removed + 1;
                end
                ctr = ctr + 1;
            end
        elseif len_to_update < batch_size
            additional = batch_size - len_to_update;
            ctr = 1; added = 0;
            [~, sorted_coinflips] = sort(coinflips, 'descend');
            while added < additional
                if ~ismember(sorted_coinflips(ctr), to_update)
                    to_update = cat(1, to_update, sorted_coinflips(ctr));
                    added = added + 1;
                end
                ctr = ctr + 1;
            end
        end
    end

    to_update = to_update(:);

end