% script to test sam pounders on either a generalized rosenbrock or cube objective.  

%% To run this test, you must ensure that the following are on your MATLAB
% path: 
% /path/to/IBCDFO/minq/m/minq5
% /path/to/IBCDFO/pounders/m
% /path/to/IBCDFO/pounders/m/general_h_funs

% With that said, this is hard-coded for MY filesystem, you must edit as 
% appropriate for however you have IBCDFO installed. 
addpath('~/IBCDFO/minq/m/minq5')
addpath('~/IBCDFO/pounders/m')
addpath('~/IBCDFO/pounders/m/general_h_funs')
addpath('../'); % don't edit this one, this is just getting sam_pounders on the path. 

%% Here, choose these strings to determine the problem you're going to run 
%test_type = 'imbalanced';
test_type = 'progressive';
%test_type = 'balanced';

macro_seed = 88;

%problem = 'rosenbrock'; 
%problem = 'cube';
%problem = 'basic_cubic';
problem = 'damped_oscillator';

if strcmp(problem, 'damped_oscillator')
    % specify time
    n = 4;
    m = 200;
    t = linspace(0, 1, m);
else
    n = 8;
    if strcmp(problem, 'rosenbrock')
        m = 2 * (n - 1); 
    else
        m = n;
    end
end

%% Less recommended to play with the stuff below. 
if ~strcmp(problem, 'damped_oscillator')
    if strcmp(test_type,'imbalanced')
        alpha = ones(1,m); alpha(m/2) = 2.0^m;
    elseif strcmp(test_type,'balanced')
        alpha = ones(1,m);
    elseif strcmp(test_type,'progressive')
        alpha = 2.0.^(1:m);
    end
else
    % define ground truth
    truth = zeros(1, n);
    truth(1) = 1.5; %A - amplitude
    truth(2) = 0.2; %gamma - damping coefficient
    truth(3) = 2.0; %omega - angular frequency
    truth(4) = 0.3; %phi - phase
    experimental_values = damped_oscillator(truth, 1:m, t);
end

figure; plot(t, experimental_values);

rng(macro_seed);

%Low = -2 * ones(1, n);
%Upp = 2 * ones(1, n); 
Low = -Inf * ones(1, n);
Upp = Inf * ones(1, n); 


if strcmp(problem,'cube')
    X0 = zeros(1,n);
elseif strcmp(problem,'rosenbrock')
    X0 = -1.0 * ones(1,n);
elseif strcmp(problem,'basic_cubic')
    X0 = ones(1, n);
elseif strcmp(problem, 'damped_oscillator')
    X0 = 0.5*ones(1, n); 
    Low = zeros(1, n);
    Upp = 2 * max(truth) * ones(1, n);
end
npmax = 2*n + 1;
g_tol = sqrt(eps);
delta_0 = 0.1;
nfs = 1;
xkin = 1;
printf = 1;
spsolver = 2;

%hfun = @(F)sum(F.^2);
%combine_models = @leastsquares; 

nf_max = 100*n;

Options = [];
Options.printf = 1;
% run pounders as a baseline
if strcmp(problem, 'rosenbrock')
    fun = @(x)generalized_rosenbrock(x, 1:m, alpha);
elseif strcmp(problem, 'cube')
    fun = @(x)generalized_cube(x, 1:m, alpha);
elseif strcmp(problem, 'basic_cubic')
    fun = @(x)basic_cubic(x, 1:m, alpha);
    Options.hfun = @(F)sum(F);
    combinemodels = @sum_combine;
    Options.combinemodels = combinemodels;
elseif strcmp(problem, 'damped_oscillator')
    fun = @(x)damped_oscillator_residual(x, 1:m, t, experimental_values);
    Options.delta_max = 0.1;
end

pounders(fun, X0, n, nf_max, g_tol, delta_0, m, Low, Upp, [], Options);
pause()

if strcmp(problem,'rosenbrock')
    fun = @(x,Set)generalized_rosenbrock(x,Set,alpha);
elseif strcmp(problem, 'cube')
    fun = @(x,Set)generalized_cube(x,Set,alpha);
elseif strcmp(problem, 'basic_cubic')
    fun = @(x, Set)basic_cubic(x, Set, alpha);
elseif strcmp(problem, 'damped_oscillator')
    fun = @(x, Set)damped_oscillator_residual(x, Set, t, experimental_values);
end

nf_max = 100*m*n;
batch_size = 1;

%expert_array = {@lipschitz_estimate_policy}; 
%expert_array = {@uniform_policy}; 

expert_array = {@lipschitz_estimate_policy, @uniform_policy}; 

%surrogate_array = rbf_surrogates(@(x)damped_oscillator_residual(x,1:m,t, experimental_values), n, m, Low, Upp);
%expert_array = {@(models, batch_size, sample_type, data)surrogate_informed_policy(models, batch_size, sample_type, data, surrogate_array, Low, Upp, spsolver)};
%expert_array = {@(models, batch_size, sample_type, data)surrogate_informed_policy(models, batch_size, sample_type, data, surrogate_array, Low, Upp, spsolver), @uniform_policy};
%expert_array = {@(models, batch_size, sample_type, data)surrogate_informed_policy(models, batch_size, sample_type, data, surrogate_array, Low, Upp, spsolver), @lipschitz_estimate_policy, @uniform_policy};

[X_inc, effort, models] = sam_pounders_paper(fun, X0, n, nf_max, g_tol, delta_0, m, Low, Upp, batch_size, expert_array, [], Options, []);

num_iters = length(models(1).critical_iters);
nfcount = zeros(m, num_iters);
for j = 1:m
    nfcount(j, 1:length(models(j).critical_iters)) = models(j).critical_iters;
    nfcount(j, 1:length(models(j).trial_iters)) = nfcount(j, 1:length(models(j).trial_iters)) + models(j).trial_iters;
    nfcount(j, 1:length(models(j).improve_iters)) = nfcount(j, 1:length(models(j).improve_iters)) + models(j).improve_iters;
    nfcount(j, 1:length(models(j).update_iters)) = nfcount(j, 1:length(models(j).update_iters)) + models(j).update_iters;
end
% the logic of the algorithm allows for iterations that don't have ANY
% evaluations. these are uninteresting in some sense, so we'll just remove
% them from this summary figure. 
nfcount(:,all(nfcount == 0))=[];
figure;
%spy(nfcount);
spyc(nfcount');