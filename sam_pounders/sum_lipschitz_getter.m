function sum = sum_lipschitz_getter(models)

    m = length(models);
    sum = 0;

    for j = 1:m
        sum = sum + models(j).Lip;
    end

end