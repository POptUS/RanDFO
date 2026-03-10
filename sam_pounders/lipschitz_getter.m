function Lip = lipschitz_getter(models, combinemodels)

    m = length(models);
    Lip = zeros(1, m);
    Cres = zeros(1, m);

    for j = 1:m
        Lip(j) = models(j).Lip;
        Cres(j) = models(j).Cres;
    end

    if strcmp(func2str(combinemodels), 'leastsquares')
        Lip = Lip .* Cres;
    end

end