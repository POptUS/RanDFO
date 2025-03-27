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
        F = [] % nf-dim array
        center_idx = 0 % idx in [1,nf] that specifies model center
        center_point = [] % (n x 1)-dim array, will always equal X(center_idx, :)
        nf = 0 % number of function evaluations (size(X, 1))
        Res = [] % nf-dim array
        Hres = [] % (n x n)-dim array containing residual Hessian
        Gres = [] % (n x 1)-dim array containing residual gradient
        Cres = [] % scalar
    end
    methods
        % constructor
        function obj = ComponentModel(id_tag, fun, X_init, F_init, xk_in, np_max, Pars, Low, Upp, delta_init, nf_max, outer_nf)
            obj.id_tag = id_tag;
            obj.fun = fun; 
            obj.X = X_init;
            obj.F = F_init;
            obj.center_idx = xk_in;
            [obj.nf, obj.n] = size(X_init);
            obj.Res = zeros(size(obj.F));
            obj.Hres = zeros(obj.n, obj.n);
            obj.Cres = obj.F(obj.center_idx);
            obj.np_max = np_max;
            obj.Pars = Pars;
            obj.Low = Low;
            obj.Upp = Upp;
            % call a function to build model on data:
            [obj.Cres, obj.Gres, obj.Hres] = build_model(obj, delta_init, nf_max, outer_nf);
        end

        function [obj, Cres, Gres, Hres, new_evals] = build_model(obj, delta, nf_max, outer_nf)

            new_evals = 0;

            % Determine the interpolation set.
            for i = 1:obj.nf
                D = obj.X(i, :) - obj.X(obj.center_idx, :);
                obj.Res(i) = obj.F(i) - obj.Cres - .5 * D * obj.Hres * D';
            end
            [Mdir, np, ~, Gres, Hresdel, ~] = ...
                formquad(obj.X(1:obj.nf, :), obj.Res(1:obj.nf, :), delta, obj.center_idx, obj.np_max, obj.Pars, 0);
            if np < obj.n  % Must obtain and evaluate bounded geometry points
                [Mdir, np] = bmpts(obj.X(obj.center_idx, :), Mdir(1:obj.n - np, :), obj.Low, obj.Upp, delta, obj.Pars(3));
                for i = 1:min(obj.n - np, nf_max - outer_nf)
                    obj.nf = obj.nf + 1;
                    obj.X(obj.nf, :) = min(obj.Upp, max(obj.Low, obj.X(obj.center_idx, :) + Mdir(i, :))); % Temp safeguard
                    obj.F(obj.nf) = obj.fun(obj.X(obj.nf, :));
                    new_evals = new_evals + 1;
                    if isnan(obj.F(obj.nf))
                        error('Nan encountered'); % need to handle this more cleanly with an exit flag, will engineer later. 
                    end
                    if printf
                        fprintf('Geometry point evaluated in component %4i \n', obj.id_tag);
                    end
                    D = Mdir(i, :);
                    obj.Res(obj.nf) = obj.F(obj.nf) - obj.Cres - .5 * D * obj.Hres * D';
                end
                [~, np, ~, Gres, Hresdel, ~] = ...
                    formquad(obj.X(1:obj.nf, :), obj.Res(1:obj.nf), delta, obj.center_idx, obj.np_max, obj.Pars, 0);
                if np < obj.n
                    error('Formquad failed to improve geometry'); % need to handle more cleanly with exit flag. 
                end
            end
        
            % 1b. Update the quadratic model
            Cres = obj.F(obj.center_idx);
            Hres = obj.Hres + Hresdel;
            obj.Cres = Cres;
            obj.Hres = Hres;
        end
    end
end