%% IMPORTANT NOTE:
% The 40 problems considered by this script are available via
% https://github.com/POptUS/YATSOp
% YATSOp should be automatically loaded in this repo as a submodule.
%
%% [1] STARS paper:
%%    https://doi.org/10.1137/22M1524072
%%    https://arxiv.org/abs/2207.06452
%% To cite the STARS paper, see Readme.md
%% at https://github.com/POptUS/RanDFO
%
% This script runs STARS in an automated way on the 40 problems considered
% in the numerical section of [1], for various types of noise, noise levels
% and subspace dimensions. It generates solutions/stats/history files
% in a 'stars_outputs' folder, which can be used to generate data profiles,
% performance profiles, trajectory plots, convergence graphs, etc.
% Users are referred to the numerical section of [1] for more details
% on the use of this script.
%
% See MCestimate.m in the 'stars_main_files' folder to understand the output files.
% Before any update related to the creation of the output files, see MCestimate.m
% and then txt_files_generator.m
%
% The noisy versions of the problems are referred to (see 'probtypes' below)
% by 'absuniform', 'reluniform2', etc. (see calfun_sample.m in the 'YATSOp'
% folder (once downloaded!) for details).
%
% See the 'Remarks' section in stars_default_options.m
% Note in particular that the algorithm uses random subspace models and as a
% consequence, the optimization of deterministic objectives may produce
% different solutions from one run to another, depending on the problem,
% unless stomads_option.FixSeed = 1.
%
% For other general experiments, see stars_applications.m
%
% See stars_default_options.m in the 'stars_main_files' folder for details
% about "stars options".
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

stars_option.yatsop = 1;    % Warning: do not modify!!
stars_option.benchmark = 0; % Warning: do not modify!!

%% Paths for algorithm main files and creation of the folder for the outputs txt files
stars_paths;  % Warning: do not comment!!

%% Sample calling syntax for the dfo and mgh functions
calldfomidfuns; % Warning: do not comment!!
stars_option.experiments = 1i; % Warning: do not modify!!

%%
global  nprobl ind_sigma

% Defining the types of problems: see calfun_sample.m in 'YATSOp' for details
probtypes = {'smooth', 'absuniform2', 'absnormal2', 'reluniform2', 'relnormal2', ...
    'absnormal', 'absuniform', 'reluniform', 'relnormal', 'abswild', 'relwild', 'nondiff'};

pvals = [1, 2, 5, 10, 20, 30, 40, 50, 100];
sigmavals = [1, 2, 3, 4, 5, 6, 7, 8];   %  See probspecs.sigma below

%% Loop on subspace dimensions
for index_p = 2 % :length(pvals)      % Indices for subspace dimension in pvals (above)

    %% Loop on the noise levels
    for ind_sigma = 3 % :length(sigmavals)   % Indices in sigmavals (above) for the standard deviation
                        % Sigma (see probspecs.sigma below)

        % One must define ind.prob, ind.sigma, ind.seed
        ind.seed = 1;
        ind.prob = 8;   % Index corresponding to the type of problem in probtypes (above)
        %                 for example, 8 corresponds to 'reluniform'
        ind.sigma = ind_sigma;        % See probspecs.sigma below
        ind.p = index_p;              % Index corresponding to the subspace dimension

        %% Define probtype and noiselevel and truncation level
        % (These will not change across the 40 problems)
        probtype = probtypes{ind.prob};
        probspecs.sigma = 10^(-sigmavals(ind.sigma));
        probspecs.trunc = 10^16;         % Chosen so that starting point unaffected

        %% Loop on the 40 problems
        for nprobl = 1:40    %  (1 to 40 problems)
            %% Initialize the rest of the problem specifications specific to problem nprobl

            % Do not modify/delete any of the probspecs below
            probspecs.nprob = Var(nprobl, 1);  % See calldfomidfuns.m for the Var array
            probspecs.n = Var(nprobl, 2);
            probspecs.m = Var(nprobl, 3);
            probspecs.hopt = Var(nprobl, 4) / 0.0001;
            probspecs.probtype = probtype;     %        This is needed for the output files

            %% Initialize the algorithm options (see stomads_default_options.m for details)
            stars_option.DisplayOutputs = 1;
            stars_option.DisplaySolution = 0;
            stars_option.FixSeed = 0;
            stars_option.HashingParam = 1;
            stars_option.HistoryFile = 1;
            stars_option.MaxFuncEval = 100000;
            stars_option.MaxRegionRadius = 5;
            stars_option.ModelLevel = 0;
            stars_option.SubspaceDim = min(probspecs.n, pvals(ind.p));
            stars_option.SampleSize = 1;
            stars_option.SolutionFile = 0;
            stars_option.StatsFile = 0;
            stars_option.SubspaceMatrix = 3;
            stars_option.UsePreviousSamples = 1;

            %%  Get starting point and problem name
            [X0, prob] = dfoxsnew(probspecs.m, probspecs.n, probspecs.nprob); % starting point
            namestr{nprobl} = prob.name;

            %% Define functions that take column vector input
            funs.myfun = @(x)calfun_sample(x, probspecs, probtype);
            funs.mysmoothfun = @(x)calfun_sample(x, probspecs, 'smooth');

            %% Run optimization  (see stars_algorithm.m in 'stars_main_files' folder)
            X = stars_algorithm(funs, X0, stars_option, probspecs);
        end
    end
end
