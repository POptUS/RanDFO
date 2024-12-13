function Q = subspace_sampler(S, sample_size)

    n = size(S, 1);

    A = randn(n, sample_size);
    
    if ~isempty(S)
        A  = A - S * (S' * A);
    end        

    [Q, ~] = qr(A);

    Q = Q(:, 1:sample_size);
end