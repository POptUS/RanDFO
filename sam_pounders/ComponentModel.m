classdef ComponentModel
    properties
        % constants
        id_tag = [] % for now, an integer indicating "which" component function this is, could be more expressive.
        fun = [] % function handle for what this component is modeling. 
        np_max = 0 % maximum number of points in model interpolation
        Pars = [] % struct containing interpolation model parameters
        Low = [] % n-dim array of lower bounds
        Upp = [] % n-dim array of upper bounds
        n = 0 % problem dimension (size(X, 2))

        % dynamic properties
        X = [] % (nf x n)-dim array
        F = zeros(1, 2) % nf-dim array
        center_idx = 0 % idx in [1,nf] that specifies model center
        center_point = [] % (n x 1)-dim array, will always equal X(center_idx, :)
        nf = 0 % number of function evaluations (size(X, 1))
        Res = [] % nf-dim array
        Hres = [] % (n x n)-dim array containing residual Hessian
        Gres = [] % (n x 1)-dim array containing residual gradient
        Cres = 0 % scalar
        Lip = 0 % Lipschitz constant estimate
        model_delta = 0 % delta current model was trained with
        valid = 0 % validity flag
        Mdir = [] % from formquad
        np = 0 % from formquad
        saved_center_idx = [] % for tentative updates
        saved_center_point = [] % for tentative updates
    end
    methods
        % constructor
        function obj = ComponentModel(id_tag, fun, X_init, F_init, xk_in, np_max, Pars, Low, Upp, delta_init, nf_max, outer_nf)
            obj.id_tag = id_tag;
            obj.fun = fun; 
            obj.X = X_init;
            obj.F(1) = F_init;
            obj.center_idx = xk_in;
            obj.center_point = obj.X; % this currently assumes a trvival Prior
            [obj.nf, obj.n] = size(X_init);
            obj.Res = zeros(size(obj.F));
            obj.Hres = zeros(obj.n, obj.n);
            obj.Cres = obj.F(obj.center_idx);
            obj.np_max = np_max;
            obj.Pars = Pars;
            obj.Low = Low;
            obj.Upp = Upp;
            % call a function to build model on data:
            obj = evaluate_initial_points(obj, delta_init);
            obj = update_model(obj, delta_init, nf_max, outer_nf);
        end

        function obj = add_new_evals(obj, y)
            obj.nf = obj.nf + 1;
            obj.X(obj.nf, :) = y;
            obj.F(obj.nf) = obj.fun(obj.X(obj.nf, :));
            if isnan(obj.F(obj.nf))
                error('Nan encountered'); % need to handle this more cleanly with an exit flag, will engineer later. 
            end
        end

        function obj = evaluate_initial_points(obj, delta_init)
            I = eye(obj.n);
            for i = 1:(obj.n)
                obj = obj.add_new_evals(obj.center_point + delta_init * I(i, :));
                obj = obj.add_new_evals(obj.center_point - delta_init * I(i, :));
            end
        end

        function [obj, new_evals] = update_model(obj, delta, nf_max, outer_nf)

            new_evals = 0;

            % Determine the interpolation set.
            for i = 1:obj.nf
                D = obj.X(i, :) - obj.X(obj.center_idx, :);
                obj.Res(i) = obj.F(i) - obj.Cres - .5 * D * obj.Hres * D';
            end
            [Mdir, np, valid, Gres, Hresdel, ~] = ...
                formquad(obj.X(1:obj.nf, :), obj.Res(1:obj.nf)', delta, obj.center_idx, obj.np_max, obj.Pars, 0);
            if np < obj.n  % Must obtain and evaluate bounded geometry points
                [Mdir, np] = bmpts(obj.X(obj.center_idx, :), Mdir(1:obj.n - np, :), obj.Low, obj.Upp, delta, obj.Pars(3));
                for i = 1:min(obj.n - np, nf_max - outer_nf)
                    y = min(obj.Upp, max(obj.Low, obj.X(obj.center_idx, :) + Mdir(i, :))); % Temp safeguard
                    obj = add_new_evals(obj, y);
                    new_evals = new_evals + 1;
                    %if printf
                        %fprintf('Geometry point evaluated in component %4i \n', obj.id_tag);
                    %end
                    D = Mdir(i, :);
                    obj.Res(obj.nf) = obj.F(obj.nf) - obj.Cres - .5 * D * obj.Hres * D';
                end
            end
            [~, np, valid, Gres, Hresdel, ~] = ...
                formquad(obj.X(1:obj.nf, :), obj.Res(1:obj.nf)', delta, obj.center_idx, obj.np_max, obj.Pars, 0);
            if np < obj.n
                error('Formquad failed to improve geometry'); % need to handle more cleanly with exit flag. 
            end
        
            % Update the quadratic model
            obj.Cres = obj.F(obj.center_idx);
            obj.Hres = obj.Hres + Hresdel;
            obj.Gres = Gres;

            % update other model properties
            obj.center_point = obj.X(obj.center_idx, :);
            obj.model_delta = delta;

            % Both of these are OK, but I like the uncommented better
            obj.Lip = max(eig(obj.Hres));
            %obj.Lip = max(obj.Lip, max(eig(obj.Hres)));

            obj.valid = valid; obj.Mdir = Mdir; obj.np = np;
        end

        function obj = update_center_point(obj, x)
            if ~all(obj.center_point == x)
                obj.center_point = x;
                obj = add_new_evals(obj, x);
                obj.center_idx = obj.nf;
            end
        end

        function obj = tentative_update_center_point(obj, x)
            obj.saved_center_idx = obj.center_idx; 
            obj.saved_center_point = obj.center_point;
            obj = update_center_point(obj, x);
        end

        function obj = make_untentative(obj)
            obj.center_idx = obj.saved_center_idx;
            obj.center_point = obj.saved_center_point;
        end

        function val = model_value_at_point(obj, y)
            displaced = y - obj.center_point;
            val = obj.Cres + displaced * obj.Gres + 0.5 * displaced * obj.Hres * displaced';
        end

        function error_estimate = get_error_estimate(obj, x_inc, tr_delta)

            % this could be reimplemented based on stuff you know about
            % your function and should consider, for instance, any hfun
            % structure. 

            % for instance, obj.Lip should be multiplied by obj.Cres if
            % hfun is leastsquares. 

            % this term is a coarse upper bound on
            % functional contributions to a difference in function values.
            center_point_dist = norm(x_inc - obj.center_point);
            if nargin == 3
                center_point_dist = center_point_dist + tr_delta;
            end
            error_estimate = obj.Lip * center_point_dist^2;
            
            % an approximation of additional contributions from model error            
            model_error_constant = sqrt(obj.n) * obj.model_delta^2 * center_point_dist / ...
                max(obj.model_delta, center_point_dist);

            error_estimate = error_estimate + obj.Lip * model_error_constant;
        end
    end
end