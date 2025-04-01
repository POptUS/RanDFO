function is_at_center = center_checker(models, x)

    m = length(models);
    is_at_center = zeros(1, m);
    for j = 1:m
        is_at_center(j) = all(x == models(j).center_point);
    end

end