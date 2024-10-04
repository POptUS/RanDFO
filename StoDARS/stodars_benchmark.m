%% IMPORTANT NOTE:
% The 53 problems considered by this script are available in
% https://github.com/POptUS/BenDFO
% which is automatically loaded in this repo as a submodule
%
%% [1] StoDARS algorithm: https://arxiv.org/abs/2403.13320
%% [2] STARS algorithm (see below): https://doi.org/10.1137/22M1524072
%% To cite the StoDARS paper [1], see Readme.md
%
%%
% This script runs stodars in an automated way on 53 problems from
% https://github.com/POptUS/BenDFO, for various types of noise and noise
% levels. It generates solutions/stats/history files in a 'stodars_outputs'
% folder, which can be used to generate data profiles, performance profiles,
% trajectory plots, convergence graphs, etc. Users can look at the
% numerical section of [2] regarding the use of this script.
%
% See stodars_MCestimate.m in the 'stodars_main_files' folder to understand the output files.
% Before any update related to the creation of the output files, see stodars_MCestimate.m
% and then stodars_txt_files_generator.m.
%
% The noisy versions of the problems are referred to (see 'probnames' below)
% by 'absuniform', 'reluniform', etc. (see calfun.m in a 'BenDFO' subfolder for details).
%
%
% See the 'Remarks and recommendations' section in stodars_default_options.m.
% Note in particular that the algorithm uses random poll directions and as a
% consequence, the optimization of deterministic objectives may produce
% different solutions from one run to another, depending on the problem,
% unless stodars_option.FixSeed = 1.
%
% For other general experiments, see stodars_applications.m.
%
% See stodars_default_options.m in the 'stodars_main_files' folder for details
% about "stodars_option".
%
%%
% Argonne National Laboratory (USA) / Lawrence Berkeley National Laboratory (USA)
%
% K. Joseph Dzahini (https://github.com/kwassi)
% and Stefan M. Wild (https://wildsm.github.io)
% October 2024

clear;                          % Be careful!!

stodars_option.yatsop = 0;    % Warning: do not modify!!
stodars_option.benchmark = 1; % Warning: do not modify!!

%% Paths for algorithm main files and creation of the folder for the outputs txt files
stodars_paths;  % Warning: do not comment!!

%% Sample calling syntax for the dfo and mgh functions
stodars_option.experiments = 1i; % Warning: do not modify!!

%%
global  nprobl dfo numprobs ind_sigma

load dfo.dat;
numprobs = size(dfo, 1);

problems_specifications;

% Defining the types of problems
probnames = {'smooth', 'absnormal', 'absuniform', 'reluniform', ...
    'relnormal', 'abswild', 'nondiff', 'relwild',  'wild3', 'noisy3'};
sigmavals = [1, 2, 3, 4, 5, 6, 7, 8];   %  See probspecs.sigma below

%% Loop on the noise levels
for ind_sigma = 3 % :length(sigmavals)     % Indices in sigmavals (above) for the standard deviation
    %                                      Sigma (see probspecs.sigma below)

    % One must define ind.prob, ind.sigma, ind.seed
    ind.seed = 1;
    ind.prob = 4;  % Index corresponding to the type of problem in probnames (above)
    %                for example, 4 corresponds to 'reluniform'
    ind.sigma = ind_sigma;        % See probspecs.sigma below

    %% Define probtype and noiselevel
    % (These will not change across the 53 problems)
    probtype = probnames{ind.prob};
    probspecs.sigma = 10^(-sigmavals(ind.sigma));

    %% Loop on the 53 problems
    for nprobl = 1:numprobs
        %% Initialize the rest of the problem specifications specific to problem nprobl
        % Do not modify/delete any of the probspecs below
        probspecs.nprob = dfo(nprobl, 1);
        probspecs.factor_power = dfo(nprobl, 4);
        probspecs.n = dfo(nprobl, 2);
        probspecs.m = dfo(nprobl, 3);
        probspecs.probtype = probtype; % This is needed for the output files

        %% Initialize the algorithm options (see stodars_default_options.m for details)
        stodars_option.DisplayOutputs = 1;
        stodars_option.DisplaySolution = 0;
        stodars_option.FixSeed = 1;
        stodars_option.HistoryFile = 1;
        stodars_option.MaxFuncEval = 60000;
        stodars_option.SampleSize = 5;
        stodars_option.SolutionFile = 1;
        stodars_option.StatsFile = 1;
        stodars_option.SubspaceDim = 2;
        stodars_option.UsePreviousSamples = 1;

        %%  Get starting point
        X0 = dfoxs(probspecs.n, probspecs.nprob, 10^probspecs.factor_power); % starting point

        %% Define a function that takes column vector input
        funs.myfun = @(x)calfun(x, probspecs, probspecs.probtype);
        funs.mysmoothfun = @(x)calfun(x, probspecs, 'smooth');

        %% Run optimization  (see stodars_algorithm.m in 'stodars_main_files' folder)
        X = stodars_algorithm(funs, X0, stodars_option, probspecs);
    end
end
