function [X_inc_array, nf_array, models] = sam_pounders_paper(Ffun, X_0, n, nf_max, g_tol, delta_0, m, Low, Upp, batch_size, expert_array, Prior, Options, Model)

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
    Options.eta_1 = sqrt(eps);     % [dbl] Parameter for accepting point, 0<eta_1<1 (.2)
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
    %[here_path, ~, ~] = fileparts(mfilename('fullpath'));
    %addpath(fullfile(here_path, 'general_h_funs'));
    hfun = @(F)sum(F.^2);
    combinemodels = @combine_leastsquares;
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
    Model.Par = zeros(1, 5);

    Model.Par(1) = sqrt(n); % [dbl] delta multiplier for checking validity
    Model.Par(2) = max(10, sqrt(n)); % [dbl] delta multiplier for all interp. points
    Model.Par(3) = 1e-5;  % [dbl] Pivot threshold for validity (1e-5)
    Model.Par(4) = .001;  % [dbl] Pivot threshold for additional points (.001)
    Model.Par(5) = false; % [bool] reverse order is false. 
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

if spsolver == 2 % Arnold Neumaier's minq5
    [here_path, ~, ~] = fileparts(mfilename('fullpath'));
    minq_path = fullfile(here_path, '..', '..', 'minq');
    addpath(fullfile(minq_path, 'm', 'minq5'));
elseif spsolver == 3 % Arnold Neumaier's minq8
    [here_path, ~, ~] = fileparts(mfilename('fullpath'));
    minq_path = fullfile(here_path, '..', '..', 'minq');
    addpath(fullfile(minq_path, 'm', 'minq8'));
end

% 0. Check inputs
[flag, X_0, np_max, ~, Low, Upp] = ...
    checkinputss(Ffun, X_0, n, Model.np_max, nf_max, g_tol, delta, nfs, m, Prior.F_init, Prior.xk_in, Low, Upp);
xk_in = 1;
if flag == -1 % Problem with the input
    X = [];
    F = [];
    hF = [];
    return
end

nf = 0; % count component function evals performed by method. 
if nfs == 0 % Need to do the first evaluation
    % Note that we're doing this evaluation FOR DISPLAY PURPOSES. 
    % We are not actually using this evaluation, which is why it does not
    % increment an nf counter. 
    % This should be fixed/removed later. 
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
        models(j) = ComponentModel(id_tag, fun, X_0, F0(j), xk_in, np_max, Model.Par, Low, Upp, delta_0, nf_max, 1, batch_size);
        % update nf in this scope:
        nf = nf + models(j).nf;        
    end
    X_inc = X_0; % explicitly store the current incumbent
    X_inc_array = [X_0; repmat(X_0, n, 1) + delta * eye(n)];
else % Have other function values around
    % populate ComponentModels
    for id_tag = 1:m
        fun = @(x)Ffun(x, id_tag);
        % Notice that this assumes all components of F were evaluated at
        % every point in X_0. Future engineering will have to worry about
        % what to do with / whether to accept partial Fvecs at points X_0. 
        models(id_tag) = ComponentModel(id_tag, fun, X_0(1:nfs, :), Prior.F_init(1:nfs, id_tag)', xk_in, np_max, Model.Par, Low, Upp, delta_0, nf_max, 1, batch_size);
        % update nf in this scope:
        %nf = nf + models(id_tag).nf;        
    end
    X_inc = X_0(xk_in, :); % explicitly store the current incumbent
    X_inc_array = X_0;
end
nf = (n + 1) * m;
nf_array = linspace(m, nf, n + 1);
success_count = n + 1;
nf_max = nf + nf_max; 

% since we just computed all the models, we effectively just did this in 
% the last step of the main loop:
to_update = (1:m)'; 
% this thing gets reset after every TR radius update: 
already_updated = to_update; 
combined_probs = ones(1, m);
[Cres, Gres, Hres] = build_average_model(models, X_inc);

% parameters for Exp4
num_experts = length(expert_array);
%if num_experts > 1
    K = nf_max / batch_size; % this is a guess of the maximum number of Exp4 rounds that can be played within budget. 
    fudge_factor = 10; 
    Exp4gamma = fudge_factor * sqrt((m * log(max(2, num_experts))) / (batch_size * K));
%else
%    Exp4gamma = 1.0;
%end

% equal initial weights on experts by default (can/should be exposed) 
weights_model = ones(num_experts, 1) / num_experts; 
weights_rho = ones(num_experts, 1) / num_experts;

% set a parameter for estimating dynamic reward scaling
EMAweight = 0.8;
EMAc = 3.0; 

% main loop
iter = 1; % this counter is for being able to reconstruct, from the models 
% class object, the history of when each component function was evaluated,
% and why. 

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
            nf = nf + 1;
            [models(j), new_evals] = update_model(models(j), delta, nf_max, nf);
            nf = nf + new_evals;
            if nf >= nf_max
                break
            end
            models(j).critical_iters(iter) = new_evals;
        end
        if nf >= nf_max
            break
        end
        valid_models = validity_checker(models);
        % update any invalid models
        for j = find(~valid_models)
            Mdir_j = models(j).Mdir; 
            [Mdir_j, np_j] = bmpts(X_inc, Mdir_j, Low, Upp, delta, Model.Par(3));
            for nf_j = 1:(n-np_j)
                nf = nf + 1;
                models(j) = add_new_evals(models(j), X_inc + Mdir_j(nf_j, :)); 
                if length(models(j).critical_iters) == iter
                    models(j).critical_iters(iter) = models(j).critical_iters(iter) + 1;
                else
                    models(j).critical_iters(iter) = 1;
                end
                if nf >= nf_max
                    break
                end
            end
            if nf >= nf_max
                break
            end
            % now update the model
            models(j) = update_center_point(models(j), X_inc);
            [models(j), new_evals] = update_model(models(j), delta, nf_max, nf);
            nf = nf + new_evals; % new_evals should actually be 0 here, but just in case!
            if nf >= nf_max
                break
            end
        end
        if nf >= nf_max
            break
        end

        % recalculate model gradient based on updated models
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
        for i = 1:n
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

        %% Evaluate Xsp
        % what does the model predict before doing any updates?
        model_prediction_inc = zeros(1, m);
        model_prediction_sp = zeros(1, m);
        for j = 1:m
            model_prediction_inc(j) = models(j).model_value_at_point(X_inc);
            model_prediction_sp(j) = models(j).model_value_at_point(Xsp);
        end

        % Generate a new sample
        additional_context = struct();
        additional_context.X_inc = X_inc;
        additional_context.Xsp = Xsp;
        additional_context.hfun = combinemodels;
    
        probs = zeros(m, num_experts); 
        for j = 1:num_experts
            probs(:, j) = expert_array{j}(models, batch_size, 'rho', additional_context);
        end
    
        combined_probs = probs * weights_rho; 
        [cps_probs, to_update] = poisson_sample_shortcut(combined_probs, batch_size);
        
        combined_probs = cps_probs';
    
        combined_probs = (1.0 - Exp4gamma) * (combined_probs / batch_size) + Exp4gamma * ones(m, 1) / m; 
    
        if any(isnan(combined_probs))
            error('nan prob')
        end

        % Now do the evals.
        valid_models = validity_checker(models); % do this first before we tentatively update model centers.

        % AMELIORATED:
        %[models, FX_inc, FXsp, new_evals] = evaluate_ameliorated_model(models, combined_probs, to_update, X_inc, Xsp, model_prediction_inc, model_prediction_sp, iter, delta, nf, nf_max);
        % AVERAGE:
        %[models, FX_inc, FXsp, new_evals] = evaluate_average_model(models, to_update, X_inc, Xsp, delta, nf_max, nf, iter);
        [models, FX_inc, FXsp, average_FX_inc, average_FXsp, new_evals] = evaluate_two_points(models, X_inc, Xsp, delta, to_update, combined_probs, iter);

        nf = nf + new_evals;

        %% update the weights
        % How bad were your model predictions?
        for j = to_update'
            reward = abs(model_prediction_inc(j) - FX_inc(j));
            reward = max(reward, abs(model_prediction_sp(j) - FXsp(j)));
            % Update the scale.
            maxe = max(eps, reward);
            if iter == 1
                reward_scale_rho = maxe;
            else
                reward_scale_rho = EMAweight * reward_scale_rho + (1.0 - EMAweight) * maxe; 
            end
            % now update weights with reward:
            scaled_reward = reward / combined_probs(j);
            for ne = 1:num_experts
                weights_rho(ne) = weights_rho(ne) * exp(Exp4gamma * probs(j, ne) * scaled_reward / (m * (EMAc * reward_scale_rho)));
            end
        end
        % And normalize the weights.
        weights_rho = weights_rho / sum(weights_rho); 

        %% compute numerator of success ratio rho
        %numerator = hfun(FXsp) - hfun(FX_inc);
        numerator = hfun(average_FXsp) - hfun(average_FX_inc);
        rho = numerator / mdec; 
    
        %% update TR center
        
        if (rho >= eta_1 || ((rho > 0) && all(valid_models)))
            avg_hF_Xsp = hfun(average_FXsp);
            if printf
                % NOTICE THAT THIS IS GIVING THE DETERMINISTIC VALUE FOR
                % THE SAKE OF EXPERIMENTATION. THE VARIABLE fval SHOULD BE
                % CHANGED TO hfun(FXsp) FOR DEPLOYMENT. 
                %fval = hfun(Ffun(Xsp, 1:m));
                fprintf('%4i  %4i  Successful iteration  %11.5e    %11.5e   %11.5e \n', nf, iter, avg_hF_Xsp, delta, ng);
            end
            X_inc = Xsp;

            % check for monotonicity in h values - if violated, increase batch size: 
            % if avg_hF_Xsp > F_inc
            %     batch_size = min(batch_size + init_batch_size, m);
            % end
            F_inc = avg_hF_Xsp;

            % output stuff:
            success_count = success_count + 1;
            X_inc_array(success_count, :) = X_inc;
            nf_array(success_count) = nf; 
            % because we tentatively updated the centers, check if we can
            % update the model center "for free" 
            for j = to_update' %1:m
                [acceptable, ~] = probe_interpolation(models(j), delta);
                if acceptable
                    [models(j), new_evals] = update_model(models(j), delta, nf_max, nf);
                    nf = nf + new_evals; % new_evals should actually be 0 here, but just in case!
                else
                    % snap the models in to_update back to their previous
                    % center
                    models(j) = make_untentative(models(j));
                    [models(j), new_evals] = update_model(models(j), delta, nf_max, nf);
                    nf = nf + new_evals; % new_evals should actually be 0 here, but just in case!
                end
                if nf >= nf_max
                    break
                end
            end
        else
            for j = to_update'
                % snap the models in to_update back to their previous
                % center
                models(j) = make_untentative(models(j));
                [models(j), new_evals] = update_model(models(j), delta, nf_max, nf);
                nf = nf + new_evals; % new_evals should actually be 0 here, but just in case!
                if nf >= nf_max
                    break
                end
            end
        end
        if nf >= nf_max
            break
        end
    
        already_updated = find(center_checker(models, X_inc));
        %close_enough = find(close_to_center_checker(models, X_inc, delta));
        %% update TR radius
        if (rho >= eta_1)  &&  (step_norm > delta_inact * delta)
            delta = min(delta * gamma_inc, delta_max);
        else
            if rho < eta_1  && all(valid_models) && length(already_updated) == m
                delta = max(delta * gamma_dec, delta_min);
                % we no longer trust that any of the models are valid in a
                % smaller ball: 
                already_updated = [];
            end
        end
    else % Don't evaluate f at Xsp
        rho = -1; % Force yourself to do a model-improving point
    end

    %% Model improvement step
    valid_models = validity_checker(models);
    if ~all(valid_models(already_updated)) && (nf < nf_max) && (rho < eta_1)
        already_updated = already_updated(:); % for annoying batch_size = 1 case
        for j = already_updated'
            [models(j), Mdir_j, ~, valid_j] = just_check_validity(models(j), delta);
            if ~valid_j
                [Mdir_j, np_j] = bmpts(X_inc, Mdir_j, Low, Upp, delta, Model.Par(3));
                if isempty(Mdir_j)
                    Mdir_j = delta * eye(n);
                    np_j = n;
                    for nf_j = 1:np_j
                        nf = nf + 1;
                        models(j) = add_new_evals(models(j), X_inc + Mdir_j(nf_j, :)); 
                    end
                else
                    % Note we can't do the same thing as pounders (greedy
                    % selection of smallest model value), because we don't know
                    % the "model value" of the combined model by changing one
                    % component. This heuristic below just picks off the first 
                    % rows of Mdir_j and seems to work just fine. Randomization 
                    % here would make reproducibility harder. Note there are 
                    % probably smarter heuristics, like the one in pounders, 
                    % but not exactly the one in pounders.
                    % IN GENERAL, this step should actually take advantage of
                    % parallel evaluations, too, and do batch_size many of them. 
                    for nf_j = 1:min(batch_size,(n-np_j)) %1
                        nf = nf + 1;
                        models(j) = add_new_evals(models(j), X_inc + Mdir_j(nf_j, :)); 
                    end
                end
               
                % now update the model
                models(j) = update_center_point(models(j), X_inc);
                [models(j), new_evals] = update_model(models(j), delta, nf_max, nf);
                nf = nf + new_evals; % new_evals should actually be 0 here, but just in case!
                if nf >= nf_max
                    break
                end
                models(j).improve_iters(iter) = min(batch_size,(n-np_j)) + new_evals;
            end
        end
    end

    %% Prepare next iteration's ameliorated model

    % Generate a new sample
    additional_context = struct();
    additional_context.X_inc = X_inc;
    additional_context.delta = delta;
    additional_context.already_updated = already_updated; 
    additional_context.hfun = combinemodels; 

    % First - is our batch size OK? 
    % need a coarse variance estimate
    %[~, var] = lipschitz_estimate_policy(models, batch_size, 'model', additional_context);
    %if var > delta^4
    %end

    probs = zeros(m, num_experts); 
    for j = 1:num_experts
        probs(:, j) = expert_array{j}(models, batch_size, 'model', additional_context);
    end

    combined_probs = probs * weights_model; 
    [cps_probs, to_update] = poisson_sample_shortcut(combined_probs, batch_size);
    combined_probs = cps_probs';

    combined_probs = (1.0 - Exp4gamma) * (combined_probs / batch_size) + Exp4gamma * ones(m, 1) / m; 

    if any(isnan(combined_probs))
        error('combined probs contains nan')
    end

     % Save copy of average model
     [Cres, Gres, Hres] = build_average_model(models, X_inc);
    
     % update the subset of models and simultaneously compute reward 
     for j = to_update'
         pre_Cres = models(j).Cres;
         pre_Gres = models(j).Gres;
         pre_Hres = models(j).Hres;
         center_point = models(j).center_point;
         displaced = X_inc - center_point;
         models(j) = update_center_point(models(j), X_inc);       
         [models(j), new_evals] = update_model(models(j), delta, nf_max, nf);
         nf = nf + new_evals;
         if nf >= nf_max
            break
         end
         models(j).update_iters(iter) = new_evals; 
         post_Cres = models(j).Cres;
         post_Gres = models(j).Gres;
         post_Hres = models(j).Hres; 
         % compute the largest change in the model in the current TR
         % see build_average_model for reminder of why this is what it is. 
         diff_model_Cres = post_Cres - pre_Cres - displaced * pre_Gres - 0.5 * displaced * pre_Hres * displaced';
         diff_model_Gres = post_Gres - pre_Gres - pre_Hres * displaced';
         diff_model_Hres = post_Hres - pre_Hres;
         % compute max_{s\inB(0,\Delta_k)} |diff_model_Cres +
         % diff_model_Gres'*s + 0.5*s'*diff_model_Hres*s|: 
        reward = max_abs_diff(diff_model_Cres, diff_model_Gres, diff_model_Hres, delta, Low, Upp, X_inc, spsolver, n);
        % NOTE: this reward also needs to be scaled (hopefully close to
        % [0,1]) 
        % Update the scale.
        maxe = max(eps, reward);
        if iter == 1
            reward_scale = maxe;
        else
            reward_scale = EMAweight * reward_scale + (1.0 - EMAweight) * maxe; 
        end
        % now update weights with reward:
        scaled_reward = reward / combined_probs(j); 
        for ne = 1:num_experts
            weights_model(ne) = weights_model(ne) * exp(Exp4gamma * probs(j, ne) * scaled_reward / (m * (EMAc * reward_scale)));
        end
     end
     if nf >= nf_max
        break
    end
     weights_model = weights_model / sum(weights_model); % normalization
     already_updated = union(already_updated, to_update'); 
    
     %% TEMPORARILY COMMENTING, SO WE USE AVERAGE MODEL HERE. 
     % build the ameliorated model
     % Notice we are assuming Cres, Gres, Hres are populated with the
     % average model data. 
     %[Cres, Gres, Hres] = build_ameliorated_model(models, X_inc, to_update, combined_probs, Cres, Gres, Hres);

    iter = iter + 1; % update the iteration counter
end % end while
end % end function