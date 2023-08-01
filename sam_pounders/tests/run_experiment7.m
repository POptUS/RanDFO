function run_experiment7(macro_seed, problem, m)

test_types = {'balanced', 'progressive', 'imbalanced'};
n = m;
sketchsize = []; 

for j = 1:3
    test_type = test_types{j};
    [ncf_vec, fvec] = sam_pounders_test(test_type, n, m, macro_seed, problem, sketchsize);
end


% min1 = min(fvec1);
% min2 = min(fvec2);
% minval = min(min1, min2) - eps;
% figure;
% plot(ncf_vec1 / ((n + 1) * m), log10(fvec1 - minval));
% hold on
% plot(ncf_vec2 / ((n + 1) * m), log10(fvec2 - minval));