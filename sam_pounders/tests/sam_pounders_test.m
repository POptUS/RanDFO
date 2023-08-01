function [ncf_vec, fvec] = sam_pounders_test(test_type, n, m, macro_seed, problem, sketchsize)

% script to test sam pounders
addpath('../')
addpath('../general_h_funs')
addpath('../../../IBCDFO/minq/m/minq5')


if strcmp(test_type,'imbalanced')
    alpha = ones(1,n); alpha(n/2) = n;
elseif strcmp(test_type,'balanced')
    alpha = ones(1,n);
elseif strcmp(test_type,'progressive')
    alpha = 1:n;
end

%Lip = [20*ones(1, n/2 - 1), 20*n, zeros(1, n/2)];

rng(macro_seed);

if strcmp(problem,'cube')
    X0 = zeros(1,n);
elseif strcmp(problem,'rosenbrock')
    X0 = -1.0 * ones(1,n);
end
npmax = 2*n + 1;
gtol = sqrt(eps);
delta = 0.1;
nfs = 1;
xkin = 1;
L = -Inf*ones(1,n);
U = Inf*ones(1,n);
printf = 1;
spsolver = 2;
combine_models = @leastsquares; 

Lip =[];
hfun = @(F)sum(F.^2);

Hf = NaN*ones(2,50*m*n);

% first do SAM
if strcmp(problem,'rosenbrock')
    fun = @(x,Set)generalized_rosenbrock(x,Set,alpha);
else
    fun = @(x,Set)generalized_cube(x,Set,alpha);
end
F0 = fun(X0, 1:m);
nfmax = 50*m*n;
[X, F, flag, xkin, ncf_vec, have_eval] = ...
sam_pounders(fun, X0, n, npmax, nfmax, gtol, delta, nfs, m, F0, xkin, L, U, sketchsize, Lip, printf, spsolver, combine_models);
nf = length(ncf_vec);
last = 0;
for j = 1:nf
    if ncf_vec(j) > 0
        Hf(1, (last + 1):min(ncf_vec(j),nfmax)) = hfun(fun(X(j, :), 1:m));
        last = ncf_vec(j);
    end
end
% now do POUNDERS
if strcmp(problem,'rosenbrock')
    fun = @(x)generalized_rosenbrock(x,1:m,alpha);
else
    fun = @(x)generalized_cube(x,1:m,alpha);
end
F0 = fun(X0);
nfmax = 50*n;
nfs = 1;
xkin = 1;
[X, F, flag, xkin] = ...
    pounders(fun, X0, n, npmax, nfmax, gtol, delta, nfs, m, F0, xkin, L, U, printf, spsolver);
ncf_vec = m * (1 : xkin);
fvec = zeros(1, xkin);
for j = 1:xkin
    Hf(2, m*(j-1)+1 : m*j) = hfun(fun(X(j, :)));
end

% save
filename = strcat('results/pounders_compare_',problem,'_',test_type,'_',num2str(macro_seed),'_',num2str(m),'.mat');
save(filename,'Hf');

end