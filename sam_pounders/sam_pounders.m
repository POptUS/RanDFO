% POUNDerS Version 0.1,    Modified 04/9/2010. Copyright 2010
% Stefan Wild and Jorge More', Argonne National Laboratory.

function [X, F, flag, xkin, ncf_vec, have_eval] = ...
    sam_pounders(fun, X0, n, npmax, nfmax, gtol, delta, nfs, m, F0, xkin, L, U, sketchsize, Lip, printf, spsolver, combinemodels, hfun)

if ~exist('hfun', 'var')
    % Use least-squares hfun by default
    addpath('../general_h_funs/');
    hfun = @(F)sum(F.^2);
    combinemodels = @leastsquares;
end
if ~exist('spsolver', 'var')
    spsolver = 2; % Use minq5 by default
end
if ~exist('printf', 'var')
    printf = 0; % Don't print by default
end
% 0. Check inputs
[flag, X0, npmax, F0, L, U] = ...
    checkinputss(fun, X0, n, npmax, nfmax, gtol, delta, nfs, m, F0, xkin, L, U);
if flag == -1 % Problem with the input
    X = [];
    F = [];
    ncf_vec = [];
    return
end

% --INTERNAL PARAMETERS [won't be changed elsewhere, defaults in ( ) ]-----
maxdelta = min(.5 * min(U - L), 1e3 * delta); % [dbl] Maximum tr radius
mindelta = min(delta * 1e-13, gtol / 10); % [dbl] Min tr radius (technically 0)
old_delta = delta * ones(1, m);
centers = ones(1, m);
gam0 = .5;      % [dbl] Parameter in (0,1) for shrinking delta  (.5)
gam1 = 2;       % [dbl] Parameter >1 for enlarging delta   (2)
eta1 = 0.05;     % [dbl] Parameter 2 for accepting point, 0<eta1<1 (.2)
Par(1) = sqrt(n); % [dbl] delta multiplier for checking validity
Par(2) = max(10, sqrt(n)); % [dbl] delta multiplier for all interp. points
Par(3) = 1e-3;  % [dbl] Pivot threshold for validity (1e-5)
Par(4) = .001;  % [dbl] Pivot threshold for additional points (.001)

% if printf
%     disp('  nf   delta    fl  np       f0           g0       ierror');
%     progstr = '%4i %9.2e %2i %3i  %11.5e %12.4e %11.3e\n'; % Line-by-line
% end
% -------------------------------------------------------------------------

% --INTERMEDIATE VARIABLES-------------------------------------------------
% D       [dbl] [1-by-n] Generic displacement vector
% G       [dbl] [n-by-1] Model gradient at X(xkin,:)
% H       [dbl] [n-by-n] Model Hessian at X(xkin,:)
% Hdel    [dbl] [n-by-n] Change to model Hessian at X(xkin,:)
% Lows    [dbl] [1-by-n] Vector of subproblem lower bounds
% Upps    [dbl] [1-by-n] Vector of subproblem upper bounds
% Mdir    [dbl] [n-by-n] Unit row directions to improve model/geometry
% Mind    [int] [npmax-by-1] Integer vector of model interpolation indices
% Xsp     [dbl] [1-by-n] Subproblem solution
% c       [dbl] Model value at X(xkin,:)
% mdec    [dbl] Change predicted by the model, m(nf)-m(xkin)
% nf      [int] Counter for the number of function evaluations
% ng      [dbl] Norm of (projection of) G
% np      [int] Number of model interpolation points
% rho     [dbl] Ratio of actual decrease to model decrease
% valid   [log] Flag saying if model is fully linear within Par(1)*delta
% -------------------------------------------------------------------------
ncf = 0; 
if nfs == 0 % Need to do the first evaluation
    X = [X0; zeros(nfmax - 1, n)]; % Stores the point locations
    F = zeros(nfmax, m); % Stores the function values
    ncf_vec = zeros(1, nfmax);
    have_eval = false(nfmax, m); % Tracks, for each point, which components we have evaluated

    nf = 1;
    F0 = fun(X(nf, :), 1:m);
    ncf = ncf + m;
    ncf_vec(nf) = ncf;
    if length(F0) ~= m
        disp('  Error: F0 does not contain the right number of residuals');
        flag = -1;
        return
    end
    F(nf, :) = F0;
    have_eval(nf, :) = true(1, m);
    if any(isnan(F(nf, :)))
        [X, F, flag] = prepare_outputs_before_return(X, F, nf, -3);
        return
    end
    if printf
        fprintf('%4i    Initial point  %11.5e\n', nf, hfun(F(nf, :)));
    end
else % Have other function values around
    X = [X0(1:nfs, :); zeros(nfmax, n)]; % Stores the point locations
    F = [F0(1:nfs, :); zeros(nfmax, m)]; % Stores the function values
    ncf_vec = zeros(1, m); 
    have_eval = [true(nfs, m); false(nfmax, m)];
    nf = nfs;
    nfmax = nfmax + nfs;
end

Cres = F(xkin, :);
Gres = zeros(n, m);
Hres = zeros(n, n, m);
ng = NaN; % Needed for early termination, e.g., if a model is never built

% For Lipschitz estimation
if isempty(Lip)
    Lip = Inf*ones(1, m); 
    %Lip = ones(1,m);
    last_valid_g = zeros(n, m);
    last_valid_center = zeros(1, m); 
    update_lip = true;
else
    update_lip = false;
end

valid = zeros(1, m);
first_success = false;
second_success = false;

while ncf < nfmax
    nf_start = nf; % for condenser
    % Step 1: Choose a subset to update
    if ~first_success || ~second_success
        models_to_update = (1:m)'; 
        probs = ones(1, m);   
    else
         [models_to_update, probs, Vm] = ... 
             choose_subset(X, Cres, xkin, centers, old_delta, delta, combinemodels, Lip, 'm', sketchsize, ng);
    end
    % Store things for the ameliorated model update
    oldGres = Gres;
    old_centers = centers; 
    oldCres = Cres; 

    % Do the corresponding function evaluations at xkin
    need_to_evaluate = intersect(models_to_update,find(~have_eval(xkin,:)));
    F(xkin, need_to_evaluate) = fun(X(xkin, :), need_to_evaluate);
    ncf = ncf + length(need_to_evaluate);
    have_eval(xkin, models_to_update) = true;
    centers(models_to_update) = xkin;
    Cres(models_to_update) = F(xkin, models_to_update);

    old_delta(models_to_update) = delta;

    % Update Gres and get Hresdel. 
    must_update = true;
    [flag, X, F, nf, ncf, Gres, Hresdel, have_eval, valid, ~, ncf_vec] = ... 
        update_model(fun, X, F, L, U, delta, xkin, nf, ncf, ... 
        Cres, Gres, Hres, models_to_update, have_eval, centers, npmax, nfmax, Par, valid, nf_start, must_update, ncf_vec);
    if ~isempty(flag)
        break
    end

    % Update Lip if necessary
    if update_lip
        [Lip, last_valid_center, last_valid_g] = update_lipschitz(Lip, X, xkin, last_valid_center, centers, last_valid_g, Gres);
    end

    % 1b. Update the quadratic model    
    [G, H] = get_ameliorated_model(X, xkin, old_centers, oldCres, Cres, oldGres, Gres, Hres, Hresdel, models_to_update, probs, combinemodels);
    Hres = Hres + Hresdel;
    ind_Lnotbinding = and(X(xkin, :) > L, G' > 0);
    ind_Unotbinding = and(X(xkin, :) < U, G' < 0);
    ng = norm(G .* (ind_Lnotbinding + ind_Unotbinding)');

    % 2. Criticality test invoked if the projected model gradient is small
    if ng < gtol
        % Check to see if the model is valid within a region of size gtol
        delta = max(gtol, max(abs(X(xkin, :))) * eps); % Safety for tiny gtols

        models_to_update = find(~have_eval(xkin, :));
        % Do the corresponding function evaluations
        F(xkin, models_to_update) = fun(X(xkin, :), models_to_update);
        ncf = ncf + length(models_to_update);
        have_eval(xkin, models_to_update) = true(1, length(models_to_update));

        centers(models_to_update) = xkin;
        Cres(models_to_update) = F(xkin, models_to_update);
        old_delta(models_to_update) = delta;

        for j = 1:m % check all the models. 
            Xj = X(have_eval(1:nf,j), :);
            Fj = F(have_eval(1:nf,j), j);
            row = 1:nf; xkinj = find(row(have_eval(:,j)) == centers(j));      
            [Mdir, ~, validj] = ...
                formquad(Xj, Fj, delta, xkinj, npmax, Par, 1);
            if ~validj % Make model valid in this small region
                [Mdir, np] = bmpts(X(xkin, :), Mdir, L, U, delta, Par(3));
                for i = 1:min(n - np, nfmax - nf)
                    nf = nf + 1;
                    X(nf, :) = min(U, max(L, X(xkin, :) + Mdir(i, :))); % Temp safeg.
                    F(nf, j) = fun(X(nf, :), j);
                    ncf = ncf + 1;
                    Xj = cat(1, Xj, X(nf, :));
                    Fj = cat(1, Fj, F(nf, j));
                    have_eval(nf, j) = true; 
                    if any(isnan(F(nf, :)))
                        [X, F, flag] = prepare_outputs_before_return(X, F, nf, -3);
                        return
                    end
                end
                if ncf >= nfmax
                    break
                end                
            end  
            % Recalculate gradient based on a MFN model
            try
                [~, ~, ~, Gresj, Hresj, ~] = ...
                    formquad(Xj, Fj, delta, xkinj, npmax, Par, 0);
                Gres(:, j) = Gresj; Hres(:, :, j) = Hresj;
            catch
                % do nothing. 
            end
        end
        [G, H] = combinemodels(Cres, Gres, Hres);
        ind_Lnotbinding = and(X(xkin, :) > L, G' > 0);
        ind_Unotbinding = and(X(xkin, :) < U, G' < 0);
        ng = norm(G .* (ind_Lnotbinding + ind_Unotbinding)');
        if ng < gtol % We trust the small gradient norm and return
            [X, F, flag] = prepare_outputs_before_return(X, F, nf, 0);
            return
        end
        models_to_update = models_to_update';
    end

    % 3. Solve the subproblem min{G'*s+.5*s'*H*s : Lows <= s <= Upps }
    Lows = max(L - X(xkin, :), -delta);
    Upps = min(U - X(xkin, :), delta);
    if spsolver == 1 % Stefan's crappy 10line solver
        [Xsp, mdec] = bqmin(H, G, Lows, Upps);
    elseif spsolver == 2 % Arnold Neumaier's minq5
        try
        H = (H + H')/2;
        catch
        pause()
        end
        [Xsp, mdec, minq_err] = minqsw(0, G, H, Lows', Upps', 0, zeros(n, 1));
        if minq_err < 0
            [X, F, flag] = prepare_outputs_before_return(X, F, nf, -4);
            return
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

    % 4. Evaluate the function at the new point (provided mdec isn't zero with an invalid model)
    if (step_norm >= 0.01 * delta || all(valid(models_to_update))) && ~(mdec == 0 && ~all(valid(models_to_update)))

        Xsp = min(U, max(L, X(xkin, :) + Xsp));  % Temp safeguard; note Xsp is not a step anymore

        % Project if we're within machine precision
        for i = 1:n % ! This will need to be cleaned up eventually
            if U(i) - Xsp(i) < eps * abs(U(i)) && U(i) > Xsp(i) && G(i) >= 0
                Xsp(i) = U(i);
                disp('eps project!');
            elseif Xsp(i) - L(i) < eps * abs(L(i)) && L(i) < Xsp(i) && G(i) >= 0
                Xsp(i) = L(i);
                disp('eps project!');
            end
        end

        if mdec == 0 && all(valid) && all(Xsp == X(xkin, :))
            [X, F, flag] = prepare_outputs_before_return(X, F, nf, -2);
            return
        end

        % 4a: Need to decide WHICH components to evaluate. 
        % store things for ameliorated estimates

        % TRYING SIMPLER THING:

        % double check that we haven't seen this trial point before
        ind = check_if_evaluated(Xsp, X, 1, nf);
        if ind == 0
            nf = nf + 1;
            X(nf, :) = Xsp; 
            ind = nf;
        end

        [f_to_update, f_probs, Vf] = ...
            choose_subset(X, Cres, ind, centers, old_delta, delta, combinemodels, Lip, 'f', sketchsize, -mdec); 
        need_to_evaluate = intersect(f_to_update, find(~have_eval(xkin,:)));
        if ~isempty(need_to_evaluate)
            F(xkin, need_to_evaluate) = fun(X(xkin, :), need_to_evaluate);
            if any(isnan(F(nf, :)))
                [X, F, flag] = prepare_outputs_before_return(X, F, nf, -3);
                return
            end
            ncf = ncf + length(need_to_evaluate);
            ncf_vec(nf) = ncf;
            have_eval(xkin, need_to_evaluate) = true;
        end

        need_to_evaluate = intersect(f_to_update, find(~have_eval(ind,:)));
        if ~isempty(need_to_evaluate)
            F(ind, need_to_evaluate) = fun(X(ind, :), need_to_evaluate);
            ncf = ncf + length(need_to_evaluate);
            ncf_vec(nf) = ncf;
            have_eval(ind, need_to_evaluate) = true;
            if any(isnan(F(ind, :)))
                [X, F, flag] = prepare_outputs_before_return(X, F, nf, -3);
                return
            end
        end

        [f0, fs] = ... 
             compute_estimates(X, xkin, ind, Cres, F(xkin, :), F(ind, :), centers, Gres, Hres, f_to_update, f_probs, f_to_update, f_probs, combinemodels);

        if mdec ~= 0
            rho = (fs - f0) / mdec;
        else % Note: this conditional only occurs when model is valid
            if fdec == 0
                [X, F, flag] = prepare_outputs_before_return(X, F, nf, -2);
                return
            else
                rho = inf * sign(fdec);
            end
        end

        % 4a. Update the center
        if (rho >= eta1)  || ((rho > 0) && all(valid))
            if ~first_success 
                first_success = true;
            elseif first_success && ~second_success
                second_success = true;
            end  
            %  Update model to reflect new center
            xkin = ind; % Change current center
            valid = zeros(1, m);
            if printf
                fprintf('%4i    Successful iteration  %11.5e   %11.5e  %11.5e \n', ncf, fs, delta, ng);
            end
            % can I update any of the models indexed by fs_to_update "for free"?
            centers_tmp = centers; Cres_tmp = Cres;
            centers_tmp(f_to_update) = xkin;
            Cres_tmp(f_to_update) = F(xkin, f_to_update);
            must_update = false;
            [flag, X, F, nf, ncf, Gres, Hresdel, have_eval, valid, freebies, ncf_vec] = ... 
                update_model(fun, X, F, L, U, delta, xkin, nf, ncf, ... 
                Cres_tmp, Gres, Hres, f_to_update, have_eval, centers_tmp, npmax, nfmax, Par, valid, nf_start, must_update, ncf_vec);
            Cres(freebies) = Cres_tmp(freebies);
            centers(freebies) = centers_tmp(freebies);
            Hres = Hres + Hresdel;
            if ~isempty(flag)
                break
            end
        else % not a success
            % can I update any of the models indexed by fc_to_update "for free"?
            centers_tmp = centers; Cres_tmp = Cres;
            centers_tmp(f_to_update) = xkin;
            Cres_tmp(f_to_update) = F(xkin, f_to_update);
            must_update = false;
            [flag, X, F, nf, ncf, Gres, Hresdel, have_eval, valid, freebies, ncf_vec] = ... 
                update_model(fun, X, F, L, U, delta, xkin, nf, ncf, ... 
                Cres_tmp, Gres, Hres, f_to_update, have_eval, centers_tmp, npmax, nfmax, Par, valid, nf_start, must_update, ncf_vec);
            Cres(freebies) = Cres_tmp(freebies);
            centers(freebies) = centers_tmp(freebies);
            Hres = Hres + Hresdel;
            if ~isempty(flag)
                break
            end
        end

        % Update Lip if necessary
        if update_lip
            [Lip, last_valid_center, last_valid_g] = update_lipschitz(Lip, X, xkin, last_valid_center, centers, last_valid_g, Gres);
        end

        % 4b. Update the trust-region radius:
        if (rho >= eta1)  &&  (step_norm > .75 * delta) %&& all(valid(models_to_update))
            delta = min(delta * gam1, maxdelta);
        elseif all(valid(models_to_update)) 
            delta = max(delta * gam0, mindelta);
        end
    else % Don't evaluate f at Xsp
        rho = -1; % Force yourself to do a model-improving point
        if printf
            disp('Warning: skipping sp soln!---------');
        end
    end

    % 5. Evaluate a model-improving point if necessary
    if ~all(valid(models_to_update))
        Res = zeros(nf, m);
        for i = 1:nf
            D = (X(i, :) - X(xkin, :));
            Res(i, :) = F(i, :) - Cres - .5 * D * reshape(D * reshape(Hres, n, m * n), n, m);
        end
    end

    for j = models_to_update'      
        if ~(valid(j)) && (ncf < nfmax) && (rho < eta1) % Implies xkin,delta unchanged
            % Need to check because model may be valid after Xsp evaluation
            Xj = X(have_eval(1:nf,j), :);
            Fj = F(have_eval(1:nf,j), j);
            Resj = Res(have_eval(1:nf, j), j);
            row = 1:nf; xkinj = find(row(have_eval(:,j)) == centers(j)); 
            [Mdir, np, validj] = ...
                formquad(Xj, Fj, delta, xkinj, npmax, Par, 1);
            if ~validj  % ! One strategy for choosing model-improving point:
                try
                % Update model (exists because delta & xkin unchanged)
                [~, ~, ~, Gresj, Hresdelj] = ...
                   formquad(Xj, Resj, delta, xkinj, npmax, Par, 0);
                Hres(:, :, j) = Hres(:, :, j) + Hresdelj;
                % Update for modelimp; Cres unchanged b/c xkin unchanged
                [G, H] = combinemodels(Cres(j), Gresj, Hres(:, :, j));
    
                % Evaluate model-improving points to pick best one
                % ! May eventually want to normalize Mdir first for infty norm
                % Plus directions
                [Mdir1, np1] = bmpts(X(xkin, :), Mdir(1:n - np, :), L, U, delta, Par(3));
                Res_tmp = zeros(n - np1, 1);
                for i = 1:n - np1
                    D = Mdir1(i, :);
                    Res_tmp(i, 1) = D * (G + .5 * H * D');
                end
                [a1, b] = min(Res_tmp(1:n - np1, 1));
                Xsp = Mdir1(b, :);
                % Minus directions
                [Mdir1, np2] = bmpts(X(xkin, :), -Mdir(1:n - np, :), L, U, delta, Par(3));
                Res_tmp = zeros(n - np2, 1);
                for i = 1:n - np2
                    D = Mdir1(i, :);
                    Res_tmp(i, 1) = D * (G + .5 * H * D');
                end
                [a2, b] = min(Res_tmp(1:n - np2, 1));
                if a2 < a1
                    Xsp = Mdir1(b, :);
                end
    
                cand = min(U, max(L, X(xkin, :) + Xsp));
                ind = check_if_evaluated(cand, X, nf_start, nf);
                if ind == 0
                    nf = nf + 1;
                    X(nf, :) = cand; % Temp safeguard
                    ind = nf;
                end

                F(ind, j) = fun(X(ind, :), j);
                ncf = ncf + 1;
                have_eval(ind, j) = true; 
                ncf_vec(nf) = ncf;

                if any(isnan(F(ind, :)))
                    [X, F, flag] = prepare_outputs_before_return(X, F, nf, -3);
                    return
                end
                catch
                    % do nothing - i don't want to stop the method just
                    % because geometry improvement failed from a numerical issue. 
                end
            end
       end
    end
end
if printf
    disp('Number of function evals exceeded');
end
flag = ng;
