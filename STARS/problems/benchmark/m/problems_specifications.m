% This file is used in stars_benchmark.m to manage the display of results.
%
% K. Joseph Dzahini (https://github.com/kwassi)
% and Stefan M. Wild (https://wildsm.github.io)
%%
% This code and its updates are available at https://github.com/POptUS/RanDFO
%%
%  Argonne National Laboratory, Mathematics and Computer Science Division
%      Kwassi Joseph Dzahini and Stefan M. Wild, October 2022
%%

global dfo numprobs

fprintf('Problem      &       n      &       m      &      Initial-Value\n');

for npl = 1:numprobs
    probspecs.nprob = dfo(npl, 1);
    probspecs.factor_power = dfo(npl, 4);
    probspecs.n = dfo(npl, 2);
    probspecs.m = dfo(npl, 3);
    y0 = dfoxs(probspecs.n, probspecs.nprob, 10^probspecs.factor_power);
    probtype = 'smooth';

    f0 = calfun(y0, probspecs, probtype);
    fprintf('%7.7g      &      %2.7g      &      %2.7g      &      %13.7g\n', ...
        npl, probspecs.n, probspecs.m, f0);
end
