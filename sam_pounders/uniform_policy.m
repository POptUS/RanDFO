function probs = uniform_policy(models, batch_size, sample_type, data)

    if strcmp(sample_type,'model')
        already_updated = data.already_updated;
    end

    m = length(models);
    if strcmp(sample_type, 'model')
        probs = zeros(m, 1);
        for j = 1:m
            if ~ismember(j, already_updated)
                % we assign 0 probs to already_updated. 
                probs(j) = 1.0 / m;
            end
        end
    elseif strcmp(sample_type, 'rho')
         probs = (1.0 / m) * ones(m, 1);
    else
        error('sample_type input to expert_policy must be a string model or rho ')
    end

sumprobs = sum(probs);
if sumprobs > 0
    % normalize
    probs = batch_size * probs;
else
    probs = (batch_size/m) * ones(m, 1);
end

end
