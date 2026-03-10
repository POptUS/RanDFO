function [L, U] = get_bounds_for_bendfo()

    addpath('~/IBCDFO/minq/m/minq5')
    addpath('~/IBCDFO/pounders/m')
    addpath('~/IBCDFO/pounders/m/general_h_funs')
    addpath('../'); % don't edit this one, this is just getting sam_pounders on the path. 
    
    addpath('~/BenDFO/m/')
    load('~/BenDFO/data/dfo.dat')
    
    num_probs = size(dfo, 1);

    L = cell(1, num_probs); U = L;

    for np = 1:num_probs
        % parameters that are constant across the tests: 
        g_tol = 1e-6; %sqrt(eps);
        delta_0 = 0.1;
        
        Options = [];
        Options.printf = 1;

        nprob = dfo(np, 1);
        n = dfo(np, 2);
        m = dfo(np, 3);
        scale_factor = dfo(np, 4);
        simplex_grads = 50;
    
        fun = @(x)bendfo_wrapper(m, n, nprob, x);
        X0 = dfoxs(n, nprob, 10^scale_factor); 
    
        % problem-dependent parameters (depend on n)
    
        Low = -Inf * ones(1, n);
        Upp = Inf * ones(1, n); 
        nf_max = simplex_grads*n;
    
        % run pounders
        [X, ~, ~] = pounders(fun, X0, n, nf_max, g_tol, delta_0, m, Low, Upp, [], Options);

        for j = 1:n
            Low(j) = min(X(:, j)) - delta_0;
            Upp(j) = max(X(:, j)) + delta_0;
        end

        L{np} = Low;
        U{np} = Upp;
    end

end