% This function returns a matrix Q such that Q-transpose is a matrix with n
% rows and p columns.
% n and p are positive integers with p <= n.
%
% K. Joseph Dzahini (https://github.com/kwassi)
% and Stefan M. Wild (https://wildsm.github.io)
%%
% This code and its updates are available at https://github.com/POptUS/RanDFO
%%
%  Argonne National Laboratory, Mathematics and Computer Science Division
%      Kwassi Joseph Dzahini and Stefan M. Wild, October 2022
%%

function Q = Gaussian_matrix(n, p)
    Q = (1 / sqrt(p)) * randn(p, n);
end
