% This file is used in stodars_benchmark.m to manage the display of results.
%
%  Argonne National Laboratory, Mathematics and Computer Science Division (USA)
%
%  Kwassi Joseph Dzahini (https://github.com/kwassi), October 2022
%%

global dfo numprobs

fprintf('Problem      &       n      &       m      &      Initial-Value\n');

for npl = 1:numprobs
    probspecs.nprob = dfo(npl, 1);
    probspecs.factor_power = dfo(npl, 4);
    probspecs.n = dfo(npl, 2);
    probspecs.m = dfo(npl, 3);
    y0 = dfoxs(probspecs.n, probspecs.nprob, 10^probspecs.factor_power);
    f0 = calfun_sample(y0, probspecs, 'smooth');
    fprintf('%7.7g      &      %2.7g      &      %2.7g      &      %13.7g\n', ...
        npl, probspecs.n, probspecs.m, f0);
end
