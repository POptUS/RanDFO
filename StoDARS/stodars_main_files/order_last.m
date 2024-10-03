% The columns of the matrix 'vectors_set' are ordered in a decreasing order
% of the cosines of their angles with respect to 'direction' (representing
% the direction from the last success in the main algorithm).
% In fact, the larger the cosine, the more the columns point in the same
% direction as 'direction'.
%
%%
%  Argonne National Laboratory (USA) / Polytechnique Montreal (Canada)
%
%  Kwassi Joseph Dzahini (https://github.com/kwassi), September 2022
%%

function Z = order_last(vectors_set, direction)
set_transpose = vectors_set';
consine_vals = zeros(length(set_transpose(:, 1)), 1);
for i = 1:length(set_transpose(:, 1))
    consine_vals(i) = cosin_vec(set_transpose(i, :), direction);
end
set_with_cosine_on_last_column = [set_transpose consine_vals];
for i = 2:length(set_with_cosine_on_last_column(:, 1))
    j = i;
    v = set_with_cosine_on_last_column(i, end);
    w = set_with_cosine_on_last_column(i, :);
    while (j > 1) && (v > set_with_cosine_on_last_column(j - 1, end))
        set_with_cosine_on_last_column(j, :) = set_with_cosine_on_last_column(j - 1, :);
        j = j - 1;
    end
    set_with_cosine_on_last_column(j, :) = w;
end
S = set_with_cosine_on_last_column';
Z = S(1:end - 1, :);
end
