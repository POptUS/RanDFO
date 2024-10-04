function z = stodars_algorithm(funs, x0, stodars_option, probspecs)
dim = length(x0);
if nargin == 3
    probspecs = struct('Dimension', 'n');
    probspecs.Dimension = dim;
end

gamma_epsilon = stodars_option.GammaEpsilon;
sampl_sze = stodars_option.SampleSize;
lw_bounds = stodars_option.LowerBounds;
if isempty(lw_bounds)
    stodars_option.LowerBounds = -Inf * ones(dim, 1);
    lw_bounds = stodars_option.LowerBounds;
end
if isrow(lw_bounds)
    lw_bounds = lw_bounds';
end
up_bounds = stodars_option.UpperBounds;
if isempty(up_bounds)
    stodars_option.UpperBounds = Inf * ones(dim, 1);
    up_bounds = stodars_option.UpperBounds;
end
if isrow(up_bounds)
    up_bounds = up_bounds';
end
prev_sampl_sze = sampl_sze;
max_iter = stodars_option.MaxNumberIters;
step_size = stodars_option.InitStepSize;
tau = stodars_option.Tau;
subsp_dim = stodars_option.SubspaceDim;
cur_sol = x0;
if stodars_option.experiments == 1i
    probspecs.Dimension = probspecs.n;
end
if isempty(probspecs.Dimension)
    probspecs.Dimension = length(cur_sol);
end
if isrow(cur_sol)
    cur_sol = cur_sol';
end
%% Warnings
stodars_warnings;
if stodars_option.warning == 1
    z = [];
    return
end
%%
success = 0;
acceptance_flag = 0;
success_flag = 0;
failure_flag = 0;

%% To make sure rng is set for reproducibility
if (stodars_option.FixSeed == 1) && (isempty(stodars_option.SeedValue))
    stodars_option.SeedValue = floor(1 + pi);
end
if stodars_option.FixSeed == 1
    rand('state', stodars_option.SeedValue);
    randn('state', stodars_option.SeedValue);
end

%%
if stodars_option.MaxFuncEval < Inf
    max_iter = 10^14;
end
global nfEval nfEvalExceeded fEval_History fEval_Stats
nfEval = 1;
nfEvalExceeded = 0;
fEval_History = [];
fEval_Stats = [];

display_cpt = 0;
if stodars_option.DisplayOutputs == 1
    if stodars_option.experiments == 1i
        fprintf('Iteration  &  Evaluations  &  Smooth-Obj  &  Current-Obj  &         Step-size\n');
    else
        fprintf('Iteration  &  Evaluations  &  Current-Obj  &         Step-size\n');
    end
end

for iteration = 1:max_iter
    %% Search
    % (Optional, no strategies proposed yet.
    % Search strategies may include, random model search (possibly in random subspaces
    % (see e.g., https://arxiv.org/abs/2207.06452)),
    % stochastic Nelder-Mead search (see https://doi.org/10.1007/s10589-018-0016-0 and
    % https://doi.org/10.1016/j.ejor.2012.02.028)
    %% Poll set
    haar_orth_mat = haar_orthog_matrix(dim);
    if stodars_option.RandomSubspaceMode == 1
        sketching_matrix = haar_orth_mat(:, 1:subsp_dim);
    else
        sketching_matrix = eye(dim);
        subsp_dim = dim;
    end
    subsp_basis = haar_orthog_matrix(subsp_dim);
    neg_sum = -sum(subsp_basis, 2);
    subsp_pos_span_set = [subsp_basis, neg_sum / norm(neg_sum)];
    fullspace_poll_dirs = sketching_matrix * subsp_pos_span_set;
    if success == 0
        poll_set = cur_sol + step_size * fullspace_poll_dirs;
    else
        poll_set = cur_sol + step_size * ...
                order_last(fullspace_poll_dirs, last_success_direction);
    end
    % Deleting poll points violating bounds constraints
    violat_vec = sum(poll_set < lw_bounds, 1) + sum(poll_set > up_bounds, 1);
    violat_index = find(violat_vec > 0);
    if isempty(violat_index) == 0
        poll_set(:, violat_index) = [];
    end
    if isempty(poll_set)
        step_size = tau * step_size;
        continue
    end
    %% Estimates computation
    % Estimate at the current solution
    funct_estm1 = stodars_MCestimate(cur_sol, sampl_sze, funs, stodars_option);
    if funct_estm1 ~= 1i
        if success_flag == 1
            cur_obj = (sampl_sze * funct_estm1 + PrevSamplSze_Succ * funct_estm2) / ...
                (sampl_sze + PrevSamplSze_Succ);
            prev_sampl_sze = sampl_sze + PrevSamplSze_Succ;
            success_flag = 0;
            PrevSamplSze_Succ = 0;
        elseif failure_flag == 1
            cur_obj = (prev_sampl_sze * cur_obj + sampl_sze * funct_estm1) / ...
                (prev_sampl_sze + sampl_sze);
            prev_sampl_sze = prev_sampl_sze + sampl_sze;
            failure_flag = 0;
        else
            cur_obj = funct_estm1;
        end
    else
        break
    end
    % Estimates at the trial points
    for j = 1:size(poll_set, 2)
        funct_estm2 = stodars_MCestimate(poll_set(:, j), sampl_sze, funs, stodars_option);
        if funct_estm2 ~= 1i
            if funct_estm2 <= cur_obj - gamma_epsilon * step_size^2
                success = 1;
                last_success_direction = fullspace_poll_dirs(:, j);   % Last direction of success
                success_point = poll_set(:, j);
                acceptance_flag = 1;
                break
            end
        else
            nfEvalExceeded = 1;
            break
        end
    end
    if nfEvalExceeded == 1
        break
    end
    %% Checking success and failure
    if acceptance_flag == 1      %  Success
        cur_sol = success_point;
        step_size = tau^(-1) * step_size;
        acceptance_flag = 0;
        if stodars_option.UsePreviousSamples == 1
            success_flag = 1;
            PrevSamplSze_Succ = sampl_sze;
        end
    else         % Failure
        step_size = tau * step_size;
        if stodars_option.UsePreviousSamples == 1
            failure_flag = 1;
        end
    end
    %% Display (or not) true function values / estimates, evaluations, iterations
    %% and step size parameters during optimization process
    display_cpt = display_cpt + 1;
    if display_cpt == 31
        if stodars_option.DisplayOutputs == 1
            if stodars_option.experiments == 1i
                fprintf('Iteration  &  Evaluations  &  Smooth-Obj  &  Current-Obj  &         Step-size\n');
            else
                fprintf('Iteration  &  Evaluations  &  Current-Obj  &         Step-size\n');
            end
        end
        display_cpt = 0;
    end
    if stodars_option.DisplayOutputs == 1
        myfun = funs.myfun;
        if stodars_option.experiments == 1i
            smooth_fun = funs.mysmoothfun;
            fprintf('%9.7g  &  %11.7g  &  %10.5g  &  %11.6g  &     %13.7g\n', ...
                    iteration, nfEval - 1, smooth_fun(cur_sol), myfun(cur_sol), step_size);
        else
            fprintf('%9.7g  &  %11.7g  &  %11.6g  &     %13.7g\n', ...
                    iteration, nfEval - 1, myfun(cur_sol), step_size);
        end
    end
end
%%
if stodars_option.experiments == 1i  % When using stodars_yatsop.m or stodars_benchmark.m
    smooth_fun = funs.mysmoothfun;
    z = [smooth_fun(cur_sol), step_size];
else  % When using stodars_applications.m
    z = [];  % Can be modified
end

%% Display or not the best solution found

if stodars_option.DisplaySolution == 1
    if stodars_option.MaxFuncEval < Inf
        Sol_message = [' Best solution found after ', num2str(nfEval - 1), ' function evaluations: '];
    else
        Sol_message = [' Best solution found after ', num2str(iteration), ' iterations: '];
    end
    disp(Sol_message);
    disp(num2str(cur_sol));
    fprintf(' \n ');
end

%%   Do not comment/delete/modify
stodars_txt_files_generator;
end
