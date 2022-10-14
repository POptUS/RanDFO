%
% K. Joseph Dzahini (https://github.com/kwassi)
% and Stefan M. Wild (https://wildsm.github.io)
%%
% This code and its updates are available at https://github.com/POptUS/RanDFO
%%
%  Argonne National Laboratory, Mathematics and Computer Science Division
%      Kwassi Joseph Dzahini and Stefan M. Wild, October 2022
%%

if length(cur_sol) ~= length(lw_bounds)
    warning('Error defining stars_option.LowerBounds and the starting point.');
    warning('Dimensions of the problem and the lower bound vector must be consistent.');
    stars_option.warning = 1;
end

if length(cur_sol) ~= length(up_bounds)
    warning('Error defining stars_option.UpperBounds and the starting point.');
    warning('Dimensions of the problem and the upper bound vector must be consistent.');
    stars_option.warning = 1;
end

if sum(cur_sol < lw_bounds) > 0
    warning('Error defining stars_option.LowerBounds or the starting point.');
    warning('Lower bound violated by the starting point.');
    stars_option.warning = 1;
end

if sum(cur_sol > up_bounds) > 0
    warning('Error defining stars_option.UpperBounds or the starting point.');
    warning('Upper bound violated by the starting point.');
    stars_option.warning = 1;
end

if sum(lw_bounds > up_bounds) > 0
    warning('Error defining stars_option.UpperBounds and stars_option.LowerBounds');
    warning('Lower bounds must be at most equal to upper bounds.');
    stars_option.warning = 1;
end

if stars_option.SubspaceDim > length(cur_sol)
    warning('Error defining stars_option.SubspaceDim.');
    warning('Subspace dimension must be at most equal to problem dimension.');
    stars_option.warning = 1;
end

if (stars_option.SubspaceMatrix == 1) && (stars_option.HashingParam > stars_option.SubspaceDim)
    warning('Error defining stars_option.HashingParam and stars_option.SubspaceDim.');
    warning('Hashing parameter cannot exceed subspace dimension.');
    stars_option.warning = 1;
end
