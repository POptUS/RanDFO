% This script runs STARS on unconstrained problems, or problems with
% bound constraints. It aims to show users how to provide problems to the
% algorithm.
%
%% STARS algorithm: (https://arxiv.org/abs/2207.06452)
%% To cite the STARS paper, see Readme.md
%
% Applications examples are provided below.
%
% See MCestimate.m in the 'stars_main_files' folder to understand the output files.
% Before any update related to the creation of the output files, see 'MCestimate.m'
% and then 'txt_files_generator.m'.
%
% For details about 'stars_option' or for other options, see 'stars_default_options.m'
% in the 'stars_main_files' folder.
%
% See the 'Remarks' section in 'stars_default_options.m'
% Note in particular that the algorithm uses random subspace models and as a
% consequence, the optimization of deterministic objectives may produce
% different solutions from one run to another, depending on the problem,
% unless stars_option.FixSeed = 1.
%
% K. Joseph Dzahini (https://github.com/kwassi)
% and Stefan M. Wild (https://wildsm.github.io)
%%
% This code and its updates are available at https://github.com/POptUS/RanDFO
%%
%  Argonne National Laboratory, Mathematics and Computer Science Division
%      Kwassi Joseph Dzahini and Stefan M. Wild, October 2022
%%

clear;                          % Be careful!!

%% Paths for algorithm main files and creation of the folder for the outputs txt files
stars_option.yatsop = 0;    % Warning: do not modify!!
stars_option.benchmark = 0; % Warning: do not modify!!
stars_paths;                  % Warning: do not comment!!
stars_option.experiments = 0; % Warning: do not comment!!

%%  Problem specifications and algorithm options
% See comments above on stars_options

stars_option.DisplayOutputs = 1;
stars_option.DisplaySolution = 1;
stars_option.FixSeed = 1;
stars_option.HashingParam = 1;
stars_option.HistoryFile = 1;           % See comments above on output files
stars_option.MaxFuncEval = 3000;
stars_option.ModelLevel = 0;
stars_option.SampleSize = 2;            % (Very small value provided here!!)
stars_option.SolutionFile = 1;
stars_option.StatsFile = 1;
stars_option.SubspaceDim = 1;
stars_option.SubspaceMatrix = 3;
stars_option.UsePreviousSamples = 1;

%%
%                          Test problems

%%                           Problem 1

% % n dimensional problem
% % Lower and upper bounds on the variables are automatically set to
% % [-Inf, ...] and [Inf, ...], respectively, since no bounds are provided
%
% % Objective function
% % noisy using a normal distribution with mean 0 and standard deviation 0.01
% funs.myfun = @(x)(norm(x)^2 - 2 * sum(x) - 27) + 0.01 * randn(1);
%
% % Starting point
% X0 = [1000, 10000, 1570, 125];   % Can be replaced by lower or higher dimension vectors

%%                           Problem 2

% % 2 dimensional problem
%
% % Objective function (deterministic part) from https://doi.org/10.1287/ijoc.1090.0319 (Section 5),
% % noisy using a uniform distribution in [-a, a]: here a = 0.01
% funs.myfun = @(x)(2 * x(1)^6 - 12.2 * x(1)^5 + 21.2 * x(1)^4 + 6.2 * x(1)^1 - 6.4 * x(1)^3 - ...
%     4.7 * x(1)^2 + x(2)^6 - 11 * x(2)^5 + 43.3 * x(2)^4 - 10 * x(2)^1 - 74.8 * ...
%     x(2)^3 + 56.9 * x(2)^2 - 4.1 * x(1) * x(2) - 0.1 * (x(1) * x(2))^2 + 0.4 * ...
%     x(1) * x(2)^2 + 0.4 * x(1)^2 * x(2)) + (0.01 + (0.01 + 0.01) * rand(1));
%
% % Starting point
% X0 = [-5, 2];                                % Wrong dimensions of X0 will throw warning messages
%                                              % (see 'stars_warnings.m' in 'stars_main_files')
%
% % Bounds on the variables
% stars_option.LowerBounds = [-Inf, 2];      % Wrong bounds will throw warning messages
% stars_option.UpperBounds = [-1, 10];

%%                           Problem 3

% Same as Problem 2, but with the objective provided by an external script (blackbox)
% in the 'problem' folder

% Adding problem location to the MATLAB path
addpath('problems/blackbox/');

% Objective function provided by a 'blackbox'
funs.myfun = @(x)bbox(x);

% Starting point
X0 = [-5, 2];                                % Wrong dimensions of X0 will throw warning messages
%                                            % (see 'stars_warnings.m' in 'stars_main_files')

% Bounds on the variables
stars_option.LowerBounds = [-Inf, 2];      % Wrong bounds will throw warning messages
stars_option.UpperBounds = [-1, 10];

%% Run optimization
X = stars_algorithm(funs, X0, stars_option);
