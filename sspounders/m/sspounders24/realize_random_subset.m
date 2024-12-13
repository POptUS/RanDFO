function Ik = realize_random_subset(pi)

    p = length(pi);

    % just a Bernoulli draw
    Ik = find(rand(p,1) <= pi);

    % however:
    if isempty(Ik)
        % must draw one at random
        Ik = randi(p,1);
    end

end