% This function computes the cosine of the angle between d and ds
% It is used in 'order_last.m' to reorder directions.
%
%%
%  Argonne National Laboratory (USA) / Polytechnique Montreal (Canada)
%
%  Kwassi Joseph Dzahini (https://github.com/kwassi), September 2022
%%

function z = cosin_vec(d, ds)
if iscolumn(d)
    d = d';
end
if iscolumn(ds)
    ds = ds';
end
z = sum(d .* ds) / (norm(d) * norm(ds));
end
