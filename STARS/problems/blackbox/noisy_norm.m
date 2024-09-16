function z = noisy_norm(x)
    z = norm(x, 'fro') + 0.1 * (2 * rand - 1);
end
