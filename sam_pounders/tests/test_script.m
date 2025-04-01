% script to test sam pounders
addpath('../general_h_funs')
addpath('../../../IBCDFO/minq/m/minq5')
addpath('../tests/')
addpath('../../../IBCDFO/pounders/m')

test_type = 'imbalanced';
macro_seed = 888;
problem = 'rosenbrock';
n = 16; 
m = n;
batch_size = 4;
dynamic_batch = false; 
% note for future struct writing: dynamic_batch is going to want a
% parameter C for the accuracy. I'm currently just setting C = n. 

if strcmp(test_type,'imbalanced')
    alpha = ones(1,n); alpha(n/2) = n;
elseif strcmp(test_type,'balanced')
    alpha = ones(1,n);
elseif strcmp(test_type,'progressive')
    alpha = 1:n;
end

rng(macro_seed);

if strcmp(problem,'cube')
    X0 = zeros(1,n);
elseif strcmp(problem,'rosenbrock')
    X0 = -1.0 * ones(1,n);
end
npmax = 2*n + 1;
g_tol = sqrt(eps);
delta_0 = 0.1;
nfs = 1;
xkin = 1;
Low = -Inf*ones(1,n);
Upp = Inf*ones(1,n);
printf = 1;
spsolver = 2;
combine_models = @leastsquares; 
nf_max = 50*m*n;

hfun = @(F)sum(F.^2);

% run pounders as a baseline
fun = @(x)generalized_rosenbrock(x, 1:m, alpha);
Options.printf = 1;
pounders(fun, X0, n, nf_max, g_tol, delta_0, m, Low, Upp, [], Options);

if strcmp(problem,'rosenbrock')
    fun = @(x,Set)generalized_rosenbrock(x,Set,alpha);
else
    fun = @(x,Set)generalized_cube(x,Set,alpha);
end

sam_pounders(fun, X0, n, nf_max, g_tol, delta_0, m, Low, Upp, batch_size, dynamic_batch);