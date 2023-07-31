function [ncf_vec, fvec] = logistic_test(test_type, n, m, seed, solver, sketchsize)

% script to test sam pounders
addpath('../')
addpath('../general_h_funs')
addpath('../../../IBCDFO/minq/m/minq5')

if strcmp(test_type,'imbalanced')
    alpha = ones(1,m); alpha(m) = m^2;
elseif strcmp(test_type,'balanced')
    alpha = ones(1,m);
elseif strcmp(test_type,'progressive')
    alpha = (1:m);%.^2;
end

rng(seed);
data = logistic_regression_data_generator(m, n, alpha);
lambda = 0.1; 
Lip = sum((data.x_train).^2)/4 + lambda;
X0 = -ones(1, n);

if strcmp(solver,'sam')
    fun = @(x,Set)log_reg(x',Set,data.x_train,data.y_train, lambda);
    F0 = fun(X0, 1:m);
    nfmax = 50*m*(n + 1);
elseif strcmp(solver,'pounders')
    fun = @(x)log_reg(x',1:m,data.x_train,data.y_train, lambda);
    F0 = fun(X0);
    nfmax = 50*(n + 1);
end

npmax = 2*n + 1;
gtol = 1e-6;
delta = 0.1;
nfs = 1;
xkin = 1;
L = -Inf*ones(1,n);
U = Inf*ones(1,n);
printf = 1;
spsolver = 2;
combine_models = @identity_combine;
hfun = @(F)sum(F);
Lip = [];

if strcmp(solver,'sam')
    [X, F, flag, xkin, ncf_vec, have_eval] = ...
        sam_pounders(fun, X0, n, npmax, nfmax, gtol, delta, nfs, m, F0, xkin, L, U, sketchsize, Lip, printf, spsolver, combine_models, hfun);
    nf = length(ncf_vec);
    fvec = zeros(1, nf);
    for j = 1:nf
        fvec(j) = hfun(fun(X(j, :), 1:m));
    end
    %figure
    %spy(have_eval(1:nf,:));
elseif strcmp(solver,'pounders')
    [X, F, flag, xkin] = ...
        pounders(fun, X0, n, npmax, nfmax, gtol, delta, nfs, m, F0, xkin, L, U, printf, spsolver, hfun, combine_models);
    ncf_vec = m * (1 : xkin);
    fvec = zeros(1, xkin);
    for j = 1:xkin
        fvec(j) = hfun(fun(X(j, :)));
    end
end

end