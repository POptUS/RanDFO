%
%% STARS algorithm (https://arxiv.org/abs/2207.06452)
%
% (Unconstrained) Stochastic derivative-free / randomized subspace method developed by
% K. Joseph Dzahini (https://github.com/kwassi)
% and Stefan M. Wild (https://wildsm.github.io)
%%
% This code and its updates are available at https://github.com/POptUS/RanDFO
%%
%  Argonne National Laboratory, Mathematics and Computer Science Division
%      Kwassi Joseph Dzahini and Stefan M. Wild, October 2022
%%

function z = stars_algorithm(funs, x0, stars_option, probspecs)
%%
if nargin == 3
    probspecs = struct('Dimension', 'n');
    probspecs.Dimension = length(x0);
end
c_param = stars_option.InflationParam;
gamma = stars_option.Gamma;
eta1 = stars_option.EtaOne;
eta2 = stars_option.EtaTwo;
TRadius = stars_option.InitRegionRadius;
TRmax = stars_option.MaxRegionRadius;
cur_sol = x0;
prob_dim = length(x0);
sub_dim = ceil(stars_option.SubspaceDim);
MaxIter = stars_option.MaxNumberIters;
JLMType = stars_option.SubspaceMatrix;
SamplSze = stars_option.SampleSize;
lw_bounds = stars_option.LowerBounds;
if isempty(lw_bounds)
    stars_option.LowerBounds = -Inf * ones(length(x0), 1);
    lw_bounds = stars_option.LowerBounds;
end
if isrow(lw_bounds)
    lw_bounds = lw_bounds';
end
up_bounds = stars_option.UpperBounds;
if isempty(up_bounds)
    stars_option.UpperBounds = Inf * ones(length(x0), 1);
    up_bounds = stars_option.UpperBounds;
end
if isrow(up_bounds)
    up_bounds = up_bounds';
end
PrevSamplSze = SamplSze;
if stars_option.experiments == 1i
    probspecs.Dimension = probspecs.n;
end
if isempty(probspecs.Dimension)
    probspecs.Dimension = length(cur_sol);
end
if isrow(cur_sol)
    cur_sol = cur_sol';
end
%% Warnings
stars_warnings;
if stars_option.warning == 1
    z = [];
    return
end

%%
hash_param = ceil(stars_option.HashingParam);
Success_flag = 0;
Failure_flag = 0;

%% To make sure rng is set for reproducibility
if (stars_option.FixSeed == 1) && (isempty(stars_option.SeedValue))
    stars_option.SeedValue = floor(1 + pi);
end
if stars_option.FixSeed == 1
    rand('state', stars_option.SeedValue);
    randn('state', stars_option.SeedValue);
end

%%
if stars_option.MaxFuncEval < Inf
    MaxIter = 10^14;
end

global nfEval nfEvalExceeded fEval_History fEval_Stats
nfEval = 1;
nfEvalExceeded = 0;
fEval_History = [];
fEval_Stats = [];

if stars_option.DisplayOutputs == 1
    if stars_option.experiments == 1i
        fprintf('Iteration  &  Evaluations  &  Smooth-Obj  &  Current-Obj  &         Radius\n');
    else
        fprintf('Iteration  &  Evaluations  &  Current-Obj  &         Radius\n');
    end
end

for iteration = 1:MaxIter
    if JLMType == 0
        subspace_matrix = Gaussian_matrix(prob_dim, sub_dim)';
    elseif JLMType == 2
        sub_dim = prob_dim;    % Q is square and full space is considered
        subspace_matrix = eye(sub_dim);
    elseif JLMType == 3
        haar_matrix = haar_orthog_matrix(prob_dim);
        subspace_matrix = sqrt(prob_dim / sub_dim) * haar_matrix(:, 1:sub_dim);
    elseif JLMType == 4
        haar_matrix = haar_orthog_matrix(prob_dim);
        subspace_matrix = haar_matrix(:, 1:sub_dim); % Without the factor sqrt(prob_dim/sub_dim)
    else
        subspace_matrix = hashing_matrix(hash_param, prob_dim, sub_dim)';
    end

    %% Approximation of the objective function value at the current solution
    % (required to estimate the gradient below)
    funct_estm1 = MCestimate(cur_sol, SamplSze, funs, stars_option);
    if funct_estm1 ~= 1i
        if Success_flag == 1
            CurObj = (SamplSze * funct_estm1 + PrevSamplSze_Succ * funct_estm3) / ...
                (SamplSze + PrevSamplSze_Succ);
            PrevSamplSze = SamplSze + PrevSamplSze_Succ;
            Success_flag = 0;
            PrevSamplSze_Succ = 0;
        elseif Failure_flag == 1
            CurObj = (PrevSamplSze * CurObj + SamplSze * funct_estm1) / (PrevSamplSze + SamplSze);
            PrevSamplSze = PrevSamplSze + SamplSze;
            Failure_flag = 0;
        else
            CurObj = funct_estm1;
        end
    else
        break
    end

    %% Model building
    % Interpolation points are the rows of poised_set
    poised_set = Algorithm_6_4(zeros(1, sub_dim), max(c_param * TRadius, 1e-100), CurObj);
    if stars_option.ModelLevel > 0
        ps = stars_option.ModelLevel;
        if sub_dim + ps > size(poised_set, 1)
            ps = size(poised_set, 1) - sub_dim;
        end
        poised_set = poised_set(1:(sub_dim + ps), :);
    end
    f_poised = CurObj;
    for idf = 2:size(poised_set, 1)
        aff_sub_point = cur_sol + subspace_matrix * poised_set(idf, :)';
        f_poised = [f_poised; MCestimate(aff_sub_point, SamplSze, funs, stars_option)];
        f_poised_end = f_poised(end);
        if f_poised_end == 1i
            nfEvalExceeded = 1;
            break
        end
    end
    if nfEvalExceeded == 1
        break
    end
    if length(f_poised) == size(poised_set, 2) + 1
        [~, GradApprox, Hessian] = MinFrobQuadMod(poised_set, f_poised);
    else
        [GradApprox, Hessian] = fitfroquad(poised_set, f_poised);
    end

    if isrow(GradApprox)
        GradApprox = GradApprox';
    end
    quad_model = @(s)(CurObj + GradApprox' * s + 0.5 * s' * Hessian * s);
    %% Solution of the Trust-Region subproblem
    NormGrad = norm(GradApprox);
    if NormGrad ~= 0
        if length(f_poised) == size(poised_set, 2) + 1
            TrialStep = -TRadius * GradApprox / NormGrad;
        else
            TrialStep = bqmin(Hessian, GradApprox, -TRadius * ones(size(poised_set, 2), 1), ...
                TRadius * ones(size(poised_set, 2), 1));
        end
        if isrow(TrialStep)
            TrialStep = TrialStep';
        end
        TrialPoint = cur_sol + subspace_matrix * TrialStep;
        %% Check whether trial point violates bound constraints or not
        if (sum(TrialPoint < lw_bounds) > 0) || (sum(TrialPoint > up_bounds) > 0) % ('Failure'!!)
            TRadius = TRadius / gamma;
            if stars_option.UsePreviousSamples == 1
                Failure_flag = 1;
            end
        else
            funct_estm3 = MCestimate(TrialPoint, SamplSze, funs, stars_option);
            if funct_estm3 == 1i
                break
            end
            %% Trust-Region ratio
            if length(f_poised) == size(poised_set, 2) + 1
                rho = (CurObj - funct_estm3) / (-GradApprox' * TrialStep);
            else
                rho = (CurObj - funct_estm3) / (quad_model(zeros(size(poised_set, 2), 1)) ...
                    - quad_model(TrialStep));
            end
            if (rho >= eta1) && (NormGrad  >= eta2 * TRadius)   % (Success)
                if stars_option.UsePreviousSamples == 1
                    Success_flag = 1;
                    PrevSamplSze_Succ = SamplSze;
                end
                cur_sol = TrialPoint;
                TRadius = min(TRmax, gamma * TRadius);
            else % (Failure)
                TRadius = TRadius / gamma;
                if stars_option.UsePreviousSamples == 1
                    Failure_flag = 1;
                end
            end
        end
    else % (Failure)
        TRadius = TRadius / gamma;
        if stars_option.UsePreviousSamples == 1
            Failure_flag = 1;
        end
    end
    %% Display (or not) true function values / estimates, evaluations, iterations
    %% and trust-region radius during optimization process
    if stars_option.DisplayOutputs == 1
        myfun = funs.myfun;
        if stars_option.experiments == 1i
            smooth_fun = funs.mysmoothfun;
            % Below, the value 'nfEval - 1' instead of 'nfEval' is due to the
            % fact that before exiting MCestimate.m, nfEval (the current number
            % of function evaluations) is incremented, i.e., nfEval = nfEval + 1
            fprintf('%9.7g  &  %11.7g  &  %10.5g  &  %11.6g  &  %13.7g\n', ...
                iteration, nfEval - 1, smooth_fun(cur_sol), myfun(cur_sol), TRadius);
        else
            fprintf('%9.7g  &  %11.7g  &  %11.6g  &  %13.7g\n', ...
                iteration, nfEval - 1, myfun(cur_sol), TRadius);
        end
    end
end

%%
if stars_option.experiments == 1i  % When using stars_yatsop.m or stars_benchmark.m
    smooth_fun = funs.mysmoothfun;
    z = [smooth_fun(cur_sol), TRadius];
else  % When using stars_applications.m
    z = [];  % Can be modified as needed
end

%% Display or not the best solution found

if stars_option.DisplaySolution == 1
    if stars_option.MaxFuncEval < Inf
        Sol_message = [' Best solution found after ', num2str(nfEval - 1), ...
            ' function evaluations: '];
    else
        Sol_message = [' Best solution found after ', num2str(iteration), ' iterations: '];
    end
    disp(Sol_message);
    disp(num2str(cur_sol));
    fprintf(' \n ');
end

%%   Do not comment/delete/modify
txt_files_generator;
end
