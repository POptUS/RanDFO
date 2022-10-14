% This function returns a matrix Q such that Q-transpose is a matrix with n
% rows and p columns.
% s denotes the number of nonzero elements out of p (that is, the number of
% nonzero elements per column of the matrix Q).
%
% K. Joseph Dzahini (https://github.com/kwassi)
% and Stefan M. Wild (https://wildsm.github.io)
%%
% This code and its updates are available at https://github.com/POptUS/RanDFO
%%
%  Argonne National Laboratory, Mathematics and Computer Science Division
%      Kwassi Joseph Dzahini and Stefan M. Wild, October 2022
%%

function Q = hashing_matrix(s, n, p)
Q = 2 * (rand(p, n) <= 0.5) - 1;  % Matrix of Rademacher (+/- 1) random variables
for j = 1:n
    B = Q(:, j);
    B(randperm(numel(B), p - s)) = 0;
    Q(:, j) = (1 / sqrt(s)) * B;
end
end
