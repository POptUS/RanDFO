% Argonne National Laboratory (USA) / Lawrence Berkeley National Laboratory (USA)
%
% K. Joseph Dzahini (https://github.com/kwassi)
% and Stefan M. Wild (https://wildsm.github.io)
% October 2024
%%

if length(cur_sol) ~= length(lw_bounds)
    warning('Error defining stodars_option.LowerBounds and the starting point.');
    warning('Dimensions of the problem and the lower bound vector must be consistent.');
    stodars_option.warning = 1;
end

if length(cur_sol) ~= length(up_bounds)
    warning('Error defining stodars_option.UpperBounds and the starting point.');
    warning('Dimensions of the problem and the upper bound vector must be consistent.');
    stodars_option.warning = 1;
end

if sum(cur_sol < lw_bounds) > 0
    warning('Error defining stodars_option.LowerBounds or the starting point.');
    warning('Lower bound violated by the starting point.');
    stodars_option.warning = 1;
end

if sum(cur_sol > up_bounds) > 0
    warning('Error defining stodars_option.UpperBounds or the starting point.');
    warning('Upper bound violated by the starting point.');
    stodars_option.warning = 1;
end

if sum(lw_bounds > up_bounds) > 0
    warning('Error defining stodars_option.UpperBounds and stodars_option.LowerBounds.');
    warning('Lower bounds must be at most equal to upper bounds.');
    stodars_option.warning = 1;
end
