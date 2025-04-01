function [X, F, hF, flag, xk_in] = sam_pounders(Ffun, X_0, n, nf_max, g_tol, delta_0, m, Low, Upp, batch_size, dynamic_batch, Prior, Options, Model)

% Check for missing arguments and initialize if necessary
if nargin < 14 || isempty(Model)
    Model = struct();
end
if nargin < 13 || isempty(Options)
    Options = struct();
end
if nargin < 12 || isempty(Prior)
    Prior = struct();
    Prior.nfs = 0;
    Prior.X_init = [];
    Prior.F_init = [];
    Prior.xk_in = 1;
end

if ~isstruct(Options)
    error("Options must be a struct");
end
if ~isstruct(Prior)
    error("Prior must be a struct");
end
if ~isstruct(Model)
    error("Model must be a struct");
end

if ~isfield(Options, 'delta_max')
    Options.delta_max = min(.5 * min(Upp - Low), 1e3 * delta_0); % [dbl] Maximum tr radius
end
if ~isfield(Options, 'delta_min')
    Options.delta_min = min(delta_0 * 1e-13, g_tol / 10); % [dbl] Min tr radius (technically 0)
end
if ~isfield(Options, 'gamma_dec')
    Options.gamma_dec = .5; % [dbl] Parameter in (0,1) for shrinking delta  (.5)
end
if ~isfield(Options, 'gamma_inc')
    Options.gamma_inc = 2;  % [dbl] Parameter (>=1) for enlarging delta   (2)
end
if ~isfield(Options, 'eta_1')
    Options.eta_1 = .05;     % [dbl] Parameter for accepting point, 0<eta_1<1 (.2)
end
if ~isfield(Options, 'delta_inact')
    Options.delta_inact = 0.75;
end
if ~isfield(Options, 'spsolver')
    Options.spsolver = 2;
end

if isfield(Options, 'hfun')
    hfun = Options.hfun;
    combinemodels = Options.combinemodels;
else
    % Use least-squares hfun by default
    [here_path, ~, ~] = fileparts(mfilename('fullpath'));
    addpath(fullfile(here_path, 'general_h_funs'));
    hfun = @(F)sum(F.^2);
    combinemodels = @leastsquares;
end
if ~isfield(Options, 'spsolver')
    Options.spsolver = 2; % Use minq5 by default
end
if ~isfield(Options, 'printf')
    Options.printf = 1; % Do print by default
end

if ~isfield(Model, 'np_max')
    Model.np_max = 2 * n + 1;
end
if ~isfield(Model, 'Par')
    Model.Par = zeros(1, 4);

    Model.Par(1) = sqrt(n); % [dbl] delta multiplier for checking validity
    Model.Par(2) = max(10, sqrt(n)); % [dbl] delta multiplier for all interp. points
    Model.Par(3) = 1e-3;  % [dbl] Pivot threshold for validity (1e-5)
    Model.Par(4) = .001;  % [dbl] Pivot threshold for additional points (.001)
end

nfs = Prior.nfs;

delta = delta_0;
spsolver = Options.spsolver;
delta_max = Options.delta_max;
delta_min = Options.delta_min;
gamma_dec = Options.gamma_dec;
gamma_inc = Options.gamma_inc;
eta_1 = Options.eta_1;
printf = Options.printf;
delta_inact = Options.delta_inact;

if     spsolver == 2 % Arnold Neumaier's minq5
    [here_path, ~, ~] = fileparts(mfilename('fullpath'));
    minq_path = fullfile(here_path, '..', '..', 'minq');
    addpath(fullfile(minq_path, 'm', 'minq5'));
elseif spsolver == 3 % Arnold Neumaier's minq8
    [here_path, ~, ~] = fileparts(mfilename('fullpath'));
    minq_path = fullfile(here_path, '..', '..', 'minq');
    addpath(fullfile(minq_path, 'm', 'minq8'));
end

% 0. Check inputs
[flag, X_0, np_max, F0, Low, Upp, xk_in] = ...
    checkinputss(Ffun, X_0, n, Model.np_max, nf_max, g_tol, delta, nfs, m, Prior.F_init, Prior.xk_in, Low, Upp);
if flag == -1 % Problem with the input
    X = [];
    F = [];
    hF = [];
    return
end

if nfs == 0 % Need to do the first evaluation
    nf = 1;
    F0 = Ffun(X_0, 1:m);
    if length(F0) ~= m
        disp('  Error: F0 does not contain the right number of residuals');
        flag = -1;
        return
    end
    if printf
        fprintf('%4i    Initial point  %11.5e\n', nf, hfun(F0));
    end
    % populate ComponentModels
    for j = 1:m
        % This structure assumes that each subset id_tag can/should be
        % computed independently. If computations were done in batches,
        % this loop would be done over batches that would be supplied as 
        % input. 
        id_tag = j;
        fun = @(x)Ffun(x, id_tag);
        models(j) = ComponentModel(id_tag, fun, X_0, F0(j), xk_in, np_max, Model.Par, Low, Upp, delta_0, nf_max, 1);
        % update nf in this scope:
        nf = nf + models(j).nf;        
    end
    X_inc = X_0; % explicitly store the current incumbent
else % Have other function values around
    nf = nfs;
    nf_max = nf_max + nfs;
    % populate ComponentModels
    for id_tag = 1:m
        fun = @(x)Ffun(x, id_tag);
        % Notice that this assumes all components of F were evaluated at
        % every point in X_0. Future engineering will have to worry about
        % what to do with / whether to accept partial Fvecs at points X_0. 
        models(j) = ComponentModel(id_tag, fun, X_0(1:nfs, :), F0(1:nfs, j), xk_in, np_max, Par, Low, Upp, delta_0, nf_max, 1);
        % update nf in this scope:
        nf = nf + models(j).nf;        
    end
    X_inc = X_0(xk_in); % explicitly store the current incumbent
end

% since we just computed all the models, we effectively did this in the
% last step of the main loop:
to_update = 1:m; 
% this thing gets reset after every TR radius change: 
already_updated = to_update; 
[Cres, Gres, Hres] = build_average_model(models, X_inc);

% main loop
while nf < nf_max

    %% Combine models
    c = hfun(Cres);
    [G, H] = combinemodels(Cres, Gres, Hres);
    ind_Lnotbinding = and(X_inc > Low, G' > 0);
    ind_Unotbinding = and(X_inc < Upp, G' < 0);
    ng = norm(G .* (ind_Lnotbinding + ind_Unotbinding)');

    %% Criticality test
    if ng < g_tol
        fprintf('Criticality step entered. \n')
        delta = max(g_tol, max(abs(X_inc)) * eps);
        correct_centers = center_checker(models, X_inc);
        % force center evaluations
        for j = find(~correct_centers)
            models(j) = update_center_point(models(j), X_inc);
            [models(j), new_evals] = update_model(models(j), delta, nf_max, nf);
            nf = nf + new_evals;
        end
        valid_models = validity_checker(models);
        % update any invalid models
        for j = find(~valid_models)
            Mdir_j = models(j).Mdir; 
            [Mdir_j, np_j] = bmpts(X_inc, Mdir_j, Low, Upp, delta, Model.Par(3));
            for nf_j = 1:(n-np_j)
                nf = nf + 1;
                models(j) = add_new_evals(models(j), X_inc + Mdir_j(nf_j, :)); 
                if nf >= nf_max
                    break
                end
            end
            % now update the model
            models(j) = update_center_point(models(j), X_inc);
            [models(j), new_evals] = update_model(models(j), delta, nf_max, nf);
            nf = nf + new_evals; % new_evals should actually be 0 here, but just in case!
            if nf >= nf_max
                break
            end
        end

        % recalculate gradient based on updated model
        [Cres, Gres, Hres] = build_average_model(models, X_inc);
        [G, H] = combinemodels(Cres, Gres, Hres);
        ind_Lnotbinding = and(X_inc > Low, G' > 0);
        ind_Unotbinding = and(X_inc < Upp, G' < 0);
        ng = norm(G .* (ind_Lnotbinding + ind_Unotbinding)');
        if ng < g_tol
            fprintf('Terminated successfully with small gradient.')
            return
        end
    end

    %% Solve TRSP
    Lows = max(Low - X_inc, -delta);
    Upps = min(Upp - X_inc, delta);
    if spsolver == 1 % Stefan's crappy 10line solver
        [Xsp, mdec] = bqmin(H, G, Lows, Upps);
    elseif spsolver == 2 % Arnold Neumaier's minq5
        H = (H + H')/2;
        [Xsp, mdec, minq_err] = minqsw(0, G, H, Lows', Upps', 0, zeros(n, 1));
        if minq_err < 0
            error('MINQ failed.')
        end
    elseif spsolver == 3 % Arnold Neumaier's minq8
        data.gam = 0;
        data.c = G;
        data.b = zeros(n, 1);
        [tmp1, tmp2] = ldl(H);
        data.D = diag(tmp2);
        data.A = tmp1';
        [Xsp, mdec] = minq8(data, Lows', Upps', zeros(n, 1), 10 * n);
    end
    Xsp = Xsp'; % Solvers currently work with column vectors
    step_norm = norm(Xsp, inf);

    valid_models = validity_checker(models);
    if (step_norm >= 0.01 * delta) || (all(valid_models) && ~(mdec == 0))

        Xsp = min(Upp, max(Low, X_inc + Xsp)); 

        % Project if we're within machine precision
        for i = 1:n % ! This will need to be cleaned up eventually
            if Upp(i) - Xsp(i) < eps * abs(Upp(i)) && Upp(i) > Xsp(i) && G(i) >= 0
                Xsp(i) = Upp(i);
                disp('eps project!');
            elseif Xsp(i) - Low(i) < eps * abs(Low(i)) && Low(i) < Xsp(i) && G(i) >= 0
                Xsp(i) = Low(i);
                disp('eps project!');
            end
        end

        valid_models = validity_checker(models);
        if mdec == 0 && all(valid_models) && all(Xsp == X_inc)
            error('No model decrease with a totally valid model!')
        end

        %% Choose a subset on which to evaluate trial Xsp
        % save the models you just updated in case you enter a model
        % improvement step at the end of this iteration
        models_to_update = to_update;

        if dynamic_batch
            sumLip = sum_lipschitz_getter(models);
            variance_estimate = Inf; this_batch_size = batch_size;
            while variance_estimate > delta^4 * n * sumLip^2
                [probs, to_update, variance_estimate] = choose_subset(models, this_batch_size, X_inc, delta, already_updated);
                this_batch_size = min(m, this_batch_size + batch_size);       
            end
        else
            [probs, to_update] = choose_subset(models, batch_size, Xsp);
            this_batch_size = batch_size;
        end

        % Save some preupdate data for the ameliorated model:
        preupdated_vals_inc = zeros(1, this_batch_size);
        preupdated_vals_sp = zeros(1, this_batch_size);
    
        % tentatively update the subset of models
        ctr = 1;
        for j = to_update'
            preupdated_vals_inc(ctr) = model_value_at_point(models(j), X_inc);
            preupdated_vals_sp(ctr) = model_value_at_point(models(j), Xsp);
            models(j) = tentative_update_center_point(models(j), Xsp);
            nf = nf + 1;
            [models(j), new_evals] = update_model(models(j), delta, nf_max, nf);
            nf = nf + new_evals;
            ctr = ctr + 1;
        end

        % evaluate ameliorated model at two points
        [FX_inc, FXsp] = evaluate_ameliorated_model(models, probs, to_update, X_inc, Xsp, preupdated_vals_inc, preupdated_vals_sp);
    
        %% compute numerator of success ratio rho
        numerator = hfun(FXsp) - hfun(FX_inc);
        rho = numerator / mdec; 
    
        %% update TR center
        valid_models = validity_checker(models);
        if rho >= eta_1 || ((rho > 0) && all(valid_models))
            if printf
                fprintf('%4i    Successful iteration  %11.5e   %11.5e  %11.5e \n', nf, c + numerator, delta, ng);
            end
            X_inc = Xsp;
            already_updated = intersect(to_update, find(valid_models));
        else
            for j = to_update'
                models(j) = make_untentative(models(j));
                % snap the models in to_update back to their previous
                % center
                [models(j), new_evals] = update_model(models(j), delta, nf_max, nf);
                nf = nf + new_evals; % new_evals should actually be 0 here, but just in case!
            end
        end
    
        %% update TR
        valid_models = validity_checker(models);
        if (rho >= eta_1)  &&  (step_norm > delta_inact * delta)
            delta = min(delta * gamma_inc, delta_max);
        elseif all(valid_models)
            delta = max(delta * gamma_dec, delta_min);
            % we no longer trust that any of the models are valid in a
            % smaller ball: 
            already_updated = [];
        end
    else % Don't evaluate f at Xsp
        rho = -1; % Force yourself to do a model-improving point
    end

    %% Model improvement step
    if ~all(valid_models) && (nf < nf_max) && (rho < eta_1)
        all_to_update = union(to_update, models_to_update);
        % for annoying batch_size = 1 case:
        all_to_update = all_to_update(:);
        for j = all_to_update'
            if ~valid_models(j)
                Mdir_j = models(j).Mdir; 
                [Mdir_j, np_j] = bmpts(X_inc, Mdir_j, Low, Upp, delta, Model.Par(3));
                % Note we can't do the same thing as pounders (greedy
                % selection of smallest model value), because we don't know
                % the "model value" of the combined model by changing one
                % component. This heuristic elow just picks off the first 
                % row of Mdir_j and seems to work just fine. Randomization 
                % here would make reproducibility harder. Note there are 
                % probably smarter heuristics, like the one in pounders, 
                % but not exactly the one in pounders.
                % IN GENERAL, this step should actually take advantage of
                % parallel evaluations, too, and do batch_size many of them. 
                for nf_j = 1:min(batch_size,(n-np_j)) %1
                    nf = nf + 1;
                    models(j) = add_new_evals(models(j), X_inc + Mdir_j(nf_j, :)); 
                end
                % now update the model
                models(j) = update_center_point(models(j), X_inc);
                [models(j), new_evals] = update_model(models(j), delta, nf_max, nf);
                nf = nf + new_evals; % new_evals should actually be 0 here, but just in case!
                % just try to fix the FIRST component in each model
                % improvement attempt. 
                break
            end
        end
    end

    %% Prepare next iteration's ameliorated model
    % Determine batch size:
    if dynamic_batch
        sumLip = sum_lipschitz_getter(models);
        variance_estimate = Inf; this_batch_size = batch_size;
        while variance_estimate > delta^4 * n * sumLip^2
            [probs, to_update, variance_estimate] = choose_subset(models, this_batch_size, X_inc, delta, already_updated);
            this_batch_size = min(m, this_batch_size + batch_size);       
        end
    else
        [probs, to_update, variance_estimate] = choose_subset(models, batch_size, X_inc, delta, already_updated);
    end

    % Save copy of average model
    [Cres, Gres, Hres] = build_average_model(models, X_inc);

    % update the subset of models
    for j = to_update'
        models(j) = update_center_point(models(j), X_inc);
        [models(j), new_evals] = update_model(models(j), delta, nf_max, nf);
        nf = nf + new_evals;
    end
    already_updated = union(already_updated, to_update');

    % build the ameliorated model
    [Cres, Gres, Hres] = build_ameliorated_model(models, X_inc, to_update, probs, Cres, Gres, Hres);
        
end % end while
end % end function