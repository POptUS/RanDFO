%% Default options for the STARS algorithm (do not change/delete/modify)
%
% K. Joseph Dzahini (https://github.com/kwassi)
% and Stefan M. Wild (https://wildsm.github.io)
%%
% This code and its updates are available at https://github.com/POptUS/RanDFO
%%
%  Argonne National Laboratory, Mathematics and Computer Science Division
%      Kwassi Joseph Dzahini and Stefan M. Wild, October 2022
%
%% [1] STARS paper:
%%    https://doi.org/10.1137/22M1524072
%%    https://arxiv.org/abs/2207.06452
%% To cite the STARS paper, see Readme.md
%% at https://github.com/POptUS/RanDFO
%
% See 'Remarks' below

stars_option.DisplayOutputs = 1;    % 1 to display #iteration, #evaluations, function values and trust-region
                                    % radius during optimization process; 0 otherwise
stars_option.DisplaySolution = 1;   % 1 to display the best solution at the end of a run; 0 otherwise
stars_option.EtaOne = 0.01;         % Algorithm parameter $\eta_1$
stars_option.EtaTwo = 0.9;          % Algorithm parameter $\eta_2$
                                    % algorithm via the trust-region radius
stars_option.FixSeed = 0;           % 1 to fix random seeds; 0 otherwise
stars_option.Gamma = 2;             % Trust-Region parameter gamma
stars_option.HashingParam = 1;      % s-Hashing parameter: should be less than subspace dimension p
stars_option.HistoryFile = 1;       % 1 to generate a .txt history file; 0 otherwise (see MCestimate.m)
stars_option.InflationParam = 1;    % Algorithm parameter $c_1 >= 1$
stars_option.InitRegionRadius = 1;  % Initial Trust-Region radius
stars_option.LowerBounds = [];      % Lower bounds on the problem variables; set to -Inf in the algorithm
                                    % if no bound is provided
stars_option.MaxFuncEval = Inf;     % Maximum number of function evaluations (handled by MCestimate.m)
stars_option.MaxNumberIters = 1000; % Maximum number of iterations (considered only when
                                    % stars_option.MaxFuncEval = Inf)
stars_option.MaxRegionRadius = 5;   % Maximum Trust-Region radius
stars_option.ModelLevel = 0;        % 0 for quadratic; 2 for linear (p+2); 3 for (p+3), etc., in subspace
stars_option.SampleSize = 5;        % Sample size for estimate computation (see MCestimate.m)
stars_option.SeedValue = [];        % Seed value when stars_option.FixSeed = 1; equals floor(1+ pi) in the
                                    % algorithm if no value is provided
stars_option.SolutionFile = 0;      % 1 to generate a .txt file containing the solution of the problem;
                                    % 0 otherwise (see [4] below)
stars_option.StatsFile = 0;         % 1 to generate a .txt stats file; 0 otherwise (see MCestimate.m)
stars_option.SubspaceDim = 1;       % subspace dimension p (from the Gaussian/Hashing matrices): default value
stars_option.SubspaceMatrix = 0;    % 0 for Gaussian, 1 for Hashing, 2 for identity matrix (STORM algorithm),
                                    % 3 for Haar-based, 4 for orthonormal columns (inspired by Haar):
                                    % matrices defining the subspaces
stars_option.UpperBounds = [];      % Upper bounds on the problem variables; set to +Inf in the algorithm
                                    % if no bound is provided
stars_option.UsePreviousSamples = 1; % 1 to use available samples in cache for estimates computation;
                                    % 0 otherwise
stars_option.warning = 0;           % Do not modify
probspecs.Dimension = [];           % Initialization! Problem dimension set to that of the starting point in
                                    % the algorithm

%% Remarks

% Please, do not uncomment the following:

%% [1]

% STORM (see stars_option.SubspaceMatrix above) is available through
% https://doi.org/10.1007/s10107-017-1141-8

%% [3]

% The algorithm uses random subspace models and as a
% consequence, the optimization of deterministic objectives may produce
% different solutions from one run to another, depending on the problem,
% unless stars_option.FixSeed = 1.

%% [4]

% When 'stars_option.HistoryFile', 'stars_option.SolutionFile' and 'stars_option.StatsFile'
% are all set to 0, the algorithm automatically sets stars_option.SolutionFile = 1 to at least
% generate the solution file of the run.
