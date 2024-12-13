function [objective, X0, n, npmax, nfmax, gtol, delta, nfs, m, F0, xkin, L, U, printf, spsolver, hfun, combinemodels] = setup_bendfo_test(probs)
addpath('~/BenDFO/m/');
addpath('~/IBCDFO/pounders/m/general_h_funs/')
addpath('~/IBCDFO/minq/m/minq5/')
addpath('~/BenDFO/data/')
addpath('../');

load dfo.dat;
numprobs = size(dfo, 1);

% Global variables needed
global BenDFO

for np = probs
    BenDFO.nprob = dfo(np, 1);
    BenDFO.n = dfo(np, 2);
    n = BenDFO.n;
    npmax = 2*n + 1;
    nfmax = 50*n;
    gtol = sqrt(eps);
    BenDFO.m = dfo(np, 3);
    m = BenDFO.m;
    BenDFO.factor_power = dfo(np, 4);

    % Obtain starting vector
    X0 = dfoxs(BenDFO.n, BenDFO.nprob, 10^BenDFO.factor_power)';
    objective = @(x)calfun_wrapper(x);

    % set remaining inputs to sspounders
    delta = 0.1;
    nfs = 1;
    F0 = objective(X0);
    xkin = 1;
    L = -Inf*ones(1,n);
    U = Inf*ones(1,n);
    printf = 1;
    spsolver = 2;
    hfun = @(F)sum(F.^2);
    combinemodels = @leastsquares;
end
end
