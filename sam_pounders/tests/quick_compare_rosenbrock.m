% quick compare rosenbrock

%test_type = 'balanced';
%test_type = 'progressive';
test_type = 'imbalanced';
n = 8;
m = 8;
sketchsize = []; 
seed = 88;

% POUNDERS
solver = 'pounders';
[ncf_vec1, fvec1] = sam_pounders_test(test_type, n, m, seed, solver);

% SAM
solver = 'sam';
[ncf_vec2, fvec2] = sam_pounders_test(test_type, n, m, seed, solver, sketchsize);

min1 = min(fvec1);
min2 = min(fvec2);
minval = min(min1, min2) - eps;
figure;
plot(ncf_vec1 / ((n + 1) * m), log10(fvec1 - minval));
hold on
plot(ncf_vec2 / ((n + 1) * m), log10(fvec2 - minval));