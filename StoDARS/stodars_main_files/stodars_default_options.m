%% Default options for the StoDARS algorithm (do not change/delete/modify)
%
%%
% Argonne National Laboratory (USA) / Lawrence Berkeley National Laboratory (USA)
%
% K. Joseph Dzahini (https://github.com/kwassi)
% and Stefan M. Wild (https://wildsm.github.io)
% October 2024
%
%%
%
% See 'Remarks and recommendations' below

stodars_option.DisplayOutputs = 1;    % 1 to display #iteration, #evaluations, function values and stepsize during optimization process, 0 otherwise
stodars_option.DisplaySolution = 1;   % 1 to display the best solution at the end of a run. 0 otherwise
stodars_option.FixSeed = 0;           % 1 to fix random seeds, 0 otherwise
stodars_option.GammaEpsilon = 1;      % See stodars paper for details regarding gamma and epsilon
stodars_option.HistoryFile = 1;       % 1 to generate a .txt history file, 0 otherwise (see stodars_MCestimate.m)
stodars_option.InitStepSize = 1;      % Initial stepsize
stodars_option.LowerBounds = [];      % Lower bounds on the problem variables; set to -Inf in the algorithm if no bound is provided
stodars_option.MaxFuncEval = Inf;     % Maximum number of function evaluations (handled by stodars_MCestimate.m)
stodars_option.MaxNumberIters = 1000; % Maximum number of iterations (considered only when stodars_option.MaxFuncEval = Inf)
stodars_option.RandomSubspaceMode = 1;% 1 for random subspace variant, 0 for the full space classical SDDS using n + 1 directions 
stodars_option.SampleSize = 5;        % Sample size for estimate computation (see stodars_MCestimate.m)
stodars_option.SeedValue = [];        % Seed value when stodars_option.FixSeed = 1; equals floor(1+ pi) in the algorithm if no value is provided
stodars_option.SolutionFile = 0;      % 1 to generate a .txt file containing the solution of the problem; 0 otherwise (see [2] below)
stodars_option.StatsFile = 0;         % 1 to generate a .txt stats file, 0 otherwise (see stodars_MCestimate.m)
stodars_option.SubspaceDim = 2;       % Value of the subspace dimension
stodars_option.Tau = 0.5;             % Parameter used to increase/decrease the stepsize parameter
stodars_option.UpperBounds = [];      % Upper bounds on the problem variables; set to +Inf in the algorithm if no bound is provided
stodars_option.UsePreviousSamples = 1;% 1 to use available samples in cache for estimates computation; 0 otherwise
stodars_option.warning = 0;           % Do not modify
probspecs.Dimension = [];             % Initialization! Problem dimension set to that of the starting point in the algorithm

%% Remarks and recommendations 

% Do not uncomment the following:

%% [1]

% The algorithm uses random orthogonal polling directions and as a
% consequence, the optimization of deterministic objectives may produce 
% different solutions from one run to another, depending on the problem,
% unless stodars_option.FixSeed = 1.

%% [2]

% When 'stodars_option.HistoryFile', 'stodars_option.SolutionFile' and 'stodars_option.StatsFile'
% are all set to 0, the algorithm automatically sets stodars_option.SolutionFile = 1 to at least
% generate the solution file of the run.
