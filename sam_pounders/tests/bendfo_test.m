% script to test sam pounders on either a generalized rosenbrock or cube objective.  

%% To run this test, you must ensure that the following are on your MATLAB
% path: 
% /path/to/IBCDFO/minq/m/minq5
% /path/to/IBCDFO/pounders/m
% /path/to/IBCDFO/pounders/m/general_h_funs

% With that said, this is hard-coded for MY filesystem, you must edit as 
% appropriate for however you have IBCDFO installed. 
addpath('~/IBCDFO/minq/m/minq5')
addpath('~/IBCDFO/pounders/m')
addpath('~/IBCDFO/pounders/m/general_h_funs')
addpath('../'); % don't edit this one, this is just getting sam_pounders on the path. 

addpath('~/BenDFO/m/')
load('~/BenDFO/data/dfo.dat')
load('bendfo_bounds.mat');

num_probs = size(dfo, 1);

% parameters that are constant across the tests: 
g_tol = 1e-6; %sqrt(eps);
delta_0 = 0.1;
nfs = 1;
xkin = 1;
printf = 1;
spsolver = 2;

Options = [];
Options.printf = 1;

hfun = @(F)sum(F.^2);

num_solvers = 4;
simplex_grads = 50;

for seed_count = 2:3
    nf_max = simplex_grads * max(dfo(:, 2).*dfo(:, 3));
    H = NaN * ones(nf_max, num_probs, num_solvers); 
    effort = H; 
    %savestr = strcat('lipschitz_test_',num2str(seed_count),'.mat');
    savestr = strcat('surrogate_test_',num2str(seed_count),'.mat');
    for np = 1:num_probs
    
        try
        nprob = dfo(np, 1);
        n = dfo(np, 2);
        m = dfo(np, 3);
        scale_factor = dfo(np, 4);
    
        fun = @(x)bendfo_wrapper(m, n, nprob, x);
        X0 = dfoxs(n, nprob, 10^scale_factor); 
    
        % problem-dependent parameters (depend on n)
    
        %Low = -Inf * ones(1, n);
        %Upp = Inf * ones(1, n); 
    
        % if using surrogate:
        Low = L{np}; Upp = U{np}; 
    
        npmax = 2*n + 1;
        nf_max = simplex_grads*n;
    
        macro_seed = 88 + seed_count - 1;
        rng(macro_seed);
    
        %% SOLVER 1
        solver_count = 1;
        if seed_count == 1
            % there is no reason to run this deterministic method more than
            % once. 
            % run pounders
            [~, ~, hF] = pounders(fun, X0, n, nf_max, g_tol, delta_0, m, Low, Upp, [], Options);
        
            % populate H and effort
            H(1:length(hF), np, solver_count) = hF; 
            effort(1:length(hF), np, solver_count) = linspace(m, m * length(hF), length(hF));
        end
    
        %% SOLVER 2 - UNIFORM ONLY
        solver_count = solver_count + 1;
        fun = @(x, Set)bendfo_wrapper(m, n, nprob, x, Set);
        nf_max = simplex_grads*m*n;
        batch_size = floor(m/2);
        if seed_count == 1 % added this option for surrogate test. 
            % sam-pounders specific:
        
            macro_seed = 88 + seed_count - 1;
            rng(macro_seed);           
            
            expert_array = {@uniform_policy}; 
            
            [X_inc_array, nf_array, models] = sam_pounders_paper(fun, X0, n, nf_max, g_tol, delta_0, m, Low, Upp, batch_size, expert_array, [], Options, []);
        
            % populate H and effort
            nX = size(X_inc_array, 1);
            for j = 1:nX
                H(j, np, solver_count) = hfun(fun(X_inc_array(j, :), 1:m));
            end
            effort(1:nX, np, solver_count) = nf_array;
        end
    
         %% SOLVER 3 - LIPSCHITZ (SURROGATE) ONLY
        solver_count = solver_count + 1;
    
        macro_seed = 88 + seed_count - 1;
        rng(macro_seed);
        
        %expert_array = {@lipschitz_estimate_policy}; 
        
        surrogate_array = rbf_surrogates(@(x)bendfo_wrapper(m, n, nprob, x, 1:m), n, m, Low, Upp);
        expert_array = {@(models, batch_size, sample_type, data)surrogate_informed_policy(models, batch_size, sample_type, data, surrogate_array, Low, Upp, spsolver)};
        
        [X_inc_array, nf_array, models] = sam_pounders_paper(fun, X0, n, nf_max, g_tol, delta_0, m, Low, Upp, batch_size, expert_array, [], Options, []);
    
        % populate H and effort
        nX = size(X_inc_array, 1);
        for j = 1:nX
            H(j, np, solver_count) = hfun(fun(X_inc_array(j, :), 1:m));
        end
        effort(1:nX, np, solver_count) = nf_array;
    
         %% SOLVER 4 - MIXED DISTRIBUTIONS
        solver_count = solver_count + 1;
    
        macro_seed = 88 + seed_count - 1;
        rng(macro_seed);
        
        %expert_array = {@lipschitz_estimate_policy, @uniform_policy}; 
        
        expert_array = {@(models, batch_size, sample_type, data)surrogate_informed_policy(models, batch_size, sample_type, data, surrogate_array, Low, Upp, spsolver), @uniform_policy};
        
        [X_inc_array, nf_array, models] = sam_pounders_paper(fun, X0, n, nf_max, g_tol, delta_0, m, Low, Upp, batch_size, expert_array, [], Options, []);
    
        % populate H and effort
        nX = size(X_inc_array, 1);
        for j = 1:nX
            H(j, np, solver_count) = hfun(fun(X_inc_array(j, :), 1:m));
        end
        effort(1:nX, np, solver_count) = nf_array;
        catch
            fprintf('This one failed. \n')
        end
        save(savestr,'H','effort','-mat');
    end
end

%% POST PROCESSING STUFF TO BE PUT BACK LATER
% num_iters = length(models(1).critical_iters);
% nfcount = zeros(m, num_iters);
% for j = 1:m
%     nfcount(j, 1:length(models(j).critical_iters)) = models(j).critical_iters;
%     nfcount(j, 1:length(models(j).trial_iters)) = nfcount(j, 1:length(models(j).trial_iters)) + models(j).trial_iters;
%     nfcount(j, 1:length(models(j).improve_iters)) = nfcount(j, 1:length(models(j).improve_iters)) + models(j).improve_iters;
%     nfcount(j, 1:length(models(j).update_iters)) = nfcount(j, 1:length(models(j).update_iters)) + models(j).update_iters;
% end
% % the logic of the algorithm allows for iterations that don't have ANY
% % evaluations. these are uninteresting in some sense, so we'll just remove
% % them from this summary figure. 
% nfcount(:,all(nfcount == 0))=[];
% figure;
% %spy(nfcount);
% spyc(nfcount');