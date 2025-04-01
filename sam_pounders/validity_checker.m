function valid_models = validity_checker(models)

    m = length(models);
    valid_models = zeros(1, m);

    for j = 1:m
        valid_models(j) = models(j).valid;
    end

end