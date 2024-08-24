%% IMPORTANT NOTE:
% The 53 problems considered by this script are available in
% https://github.com/POptUS/BenDFO
% which is automatically loaded in this repo as a submodule.
%
%% [1] STARS paper:
%%    https://doi.org/10.1137/22M1524072
%%    https://arxiv.org/abs/2207.06452
%% To cite the STARS paper, see Readme.md
%% at https://github.com/POptUS/RanDFO
%
% This script runs STARS in an automated way on 53 problems from
% https://github.com/POptUS/BenDFO, for various types of noise, noise levels
% and subspace dimensions. It generates solutions/stats/history files
% in a 'stars_outputs' folder, which can be used to generate data profiles,
% performance profiles, trajectory plots, convergence graphs, etc.
% Users can be inspired by the numerical section of [1] regarding the use
% of this script.
%
% See MCestimate.m in the 'stars_main_files' folder to understand the output files.
% Before any update related to the creation of the output files, see 'MCestimate.m'
% and then 'txt_files_generator.m'.
%
% The noisy versions of the problems are referred to (see 'probnames' below)
% by 'absuniform', 'reluniform', etc. (see 'calfun.m' in a 'BenDFO'
% subfolder (once downloaded!) for details).
%
% See the 'Remarks' section in 'stars_default_options.m'.
% Note in particular that the algorithm uses random subspace models and as a
% consequence, the optimization of deterministic objectives may produce
% different solutions from one run to another, depending on the problem,
% unless stars_option.FixSeed = 1.
%
% For other general experiments, see 'stars_applications.m'.
%
% See 'stars_default_options.m' in the 'stars_main_files' folder for details
% about "stars_option".
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

stars_option.yatsop = 0;    % Warning: do not modify!!
stars_option.benchmark = 1; % Warning: do not modify!!

%% Paths for algorithm main files and creation of the folder for the outputs txt files
stars_paths;  % Warning: do not comment!!

%% Sample calling syntax for the dfo and mgh functions

stars_option.experiments = 1i; % Warning: do not modify!!

%%

global  nprobl dfo numprobs ind_sigma

load dfo.dat;
numprobs = size(dfo, 1);

problems_specifications;

% Defining the types of problems
probnames = {'smooth', 'absnormal', 'absuniform', 'reluniform', ...
    'relnormal', 'abswild', 'nondiff', 'relwild',  'wild3', 'noisy3'};

pvals = [1, 2, 5, 10, 15, 20, 25, 31];

sigmavals = [1, 2, 3, 4, 5, 6, 7, 8];   %  See probspecs.sigma below

%% Loop on subspace dimensions
for index_p = 2 % :length(pvals)      % Indices for subspace dimension in pvals (above)

    %% Loop on the noise levels
    for ind_sigma = 3 % :length(sigmavals)   % Indices in sigmavals (above) for the standard deviation
        %                                    Sigma (see probspecs.sigma below)

        % One must define ind.prob, ind.sigma, ind.seed
        ind.seed = 1;
        ind.prob = 4;  % Index corresponding to the type of problem in probnames (above):
        %                for example, 4 corresponds to 'reluniform'
        ind.sigma = ind_sigma;        % See probspecs.sigma below
        ind.p = index_p;              % Index corresponding to the subspace dimension

        %% Define probtype and noiselevel
        %  (These will not change across the 53 problems)
        probtype = probnames{ind.prob};
        probspecs.sigma = 10^(-sigmavals(ind.sigma));

        %% Loop on the 53 problems
        for nprobl = 1:numprobs
            %% Initialize the rest of the problem specifications specific to problem nprobl
            %  Do not modify/delete any of the probspecs below
            probspecs.nprob = dfo(nprobl, 1);
            probspecs.factor_power = dfo(nprobl, 4);
            probspecs.n = dfo(nprobl, 2);
            probspecs.m = dfo(nprobl, 3);
            probspecs.hopt = stars_option.FiniDiffParam;
            probspecs.probtype = probtype; % This is needed for the output files

            %% Initialize the algorithm options (see stars_default_options.m for details)
            stars_option.DisplayOutputs = 1;
            stars_option.DisplaySolution = 0;
            stars_option.FixSeed = 1;
            stars_option.HashingParam = 1;
            stars_option.HistoryFile = 1;
            stars_option.MaxFuncEval = 6000;
            stars_option.MaxRegionRadius = 5;
            stars_option.SampleSize = 5;
            stars_option.SolutionFile = 0;
            stars_option.StatsFile = 0;
            stars_option.SubspaceMatrix = 1;
            stars_option.SubspaceDim = min(probspecs.n, pvals(ind.p));
            stars_option.UsePreviousSamples = 1;

            %%  Get starting point
            X0 = dfoxs(probspecs.n, probspecs.nprob, 10^probspecs.factor_power); % starting point

            %% Define functions that take column vector input
            funs.myfun = @(x)calfun(x, probspecs, probspecs.probtype);
            funs.mysmoothfun = @(x)calfun(x, probspecs, 'smooth');

            %% Run optimization  (see stars_algorithm.m in 'stars_main_files' folder)
            X = stars_algorithm(funs, X0, stars_option, probspecs);
        end
    end
end
