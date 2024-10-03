%% IMPORTANT NOTE:
% The 40 problems considered by this script are available via
% https://github.com/POptUS/YATSOp
% YATSOp should be automatically loaded in this repo as a submodule.
% See 'Readme.txt' in the YATSOp subfolder located in the 'problems'
% folder for other details.
%
%% [1] StoDARS algorithm: https://arxiv.org/abs/2403.13320
%% [2] STARS algorithm (see below): https://doi.org/10.1137/22M1524072
%% To cite the StoDARS paper [1], see Readme.md
%
%%
% This script runs stodars in an automated way on the 40 problems considered
% in the numerical section of [2], for various types of noise and various
% noise levels. It generates solutions files, stats files and history files
% in a 'stodars_outputs' folder, which can be used to generate data profiles,
% performance profiles, trajectory plots, convergence graphs, etc.
% Users can refer to the numerical section of the above manuscript
% for more details on the use of this script.
%
% See stodars_MCestimate.m in the 'stodars_main_files' folder to understand the output files.
% Before any update related to the creation of the output files, see 'stodars_MCestimate.m'
% and then 'stodars_txt_files_generator.m'.
%
% The noisy versions of the problems are referred to (see "probtypes" below)
% by 'absuniform', 'reluniform2', etc. (see 'calfun_sample.m' in the 'YATSOp'
% folder for details).
%
% See the 'Remarks and recommendations' section in 'stodars_default_options.m'.
% Note in particular that the algorithm uses random orthogonal polling
% directions and as a consequence, the optimization of deterministic
% objectives may produce different solutions from one run to another,
% depending on the problem, unless stodars_option.FixSeed = 1.
%
% For other general experiments, see 'stodars_applications.m'.
%
% See 'stodars_default_options.m' in the 'stodars_main_files' folder for details
% about "stodars options".
%
%%
% Argonne National Laboratory (USA) / Lawrence Berkeley National Laboratory (USA)
%
% K. Joseph Dzahini (https://github.com/kwassi)
% and Stefan M. Wild (https://wildsm.github.io)
% October 2024

clear;                          % Be careful!!

stodars_option.yatsop = 1;    % Warning: do not modify!!
stodars_option.benchmark = 0; % Warning: do not modify!!

%% Paths for algorithm main files and creation of the folder for the outputs txt files
stodars_paths;  % Warning: do not comment!!

%% Sample calling syntax for the dfo and mgh functions
calldfomidfuns; % Warning: do not comment!!
stodars_option.experiments = 1i; % Warning: do not modify!!

%%
global  nprobl ind_sigma

% Defining the types of problems: see calfun_sample.m in 'YATSOp' for details
probtypes = {'smooth', 'absuniform2', 'absnormal2', 'reluniform2', 'relnormal2', ...
    'absnormal', 'absuniform', 'reluniform', 'relnormal', 'abswild', 'relwild', 'nondiff'};

sigmavals = [1, 2, 3, 4, 5, 6, 7, 8];   %  See probspecs.sigma below

%% Loop on the noise levels
for type_val = 4:5
    for repl = 1:20
        for ind_sigma = 3 %:length(sigmavals)     % Indices in sigmavals (above) for the standard deviation
            %                                      Sigma (see probspecs.sigma below)
            % One must define ind.prob, ind.sigma, ind.seed
            ind.seed = 1;
            ind.prob = type_val;  % Index corresponding to the type of problem in probtypes (above)
            %                for example, 8 corresponds to 'reluniform'
            ind.sigma = ind_sigma;        % See probspecs.sigma below

            %% define probtype and noiselevel and truncation level
            % (these will not change across the 40 problems)
            probtype = probtypes{ind.prob};
            probspecs.sigma = 10^(-sigmavals(ind.sigma));
            probspecs.trunc = 10^16; % Chosen so that starting point unaffected

            %% Loop on the 40 problems

            for nprobl = 1:40   %  (1 to 40 problems)
                %% Initialize the rest of the problem specifications specific to problem nprobl
                % Do not modify/delete any of the probspecs below
                probspecs.nprob = Var(nprobl, 1);  % See calldfomidfuns.m for the Var array
                probspecs.n = Var(nprobl, 2);
                probspecs.m = Var(nprobl, 3);
                probspecs.probtype = probtype; % This is needed for the output files

                %% Initializing the algorithm options (see stodars_default_options.m for details)
                stodars_option.DisplayOutputs = 1;
                stodars_option.DisplaySolution = 0;
                stodars_option.FixSeed = 1;
                stodars_option.HistoryFile = 1;
                stodars_option.MaxFuncEval = 1500 * (probspecs.n + 1);
                stodars_option.RandomSubspaceMode = 1;
                stodars_option.SampleSize = 1;
                stodars_option.SeedValue = floor(1 + pi^repl);
                stodars_option.SolutionFile = 1;
                stodars_option.StatsFile = 1;
                stodars_option.UsePreviousSamples = 1;

                %%  Get starting point and problem name
                [X0, prob] = dfoxsnew(probspecs.m, probspecs.n, probspecs.nprob); % starting point
                namestr{nprobl} = prob.name;

                %% Define a function that takes column vector input
                funs.myfun = @(x)calfun_sample(x, probspecs, probtype);
                funs.mysmoothfun = @(x)calfun_sample(x, probspecs, 'smooth');

                %% Run optimization  (see stodars_algorithm.m in 'stodars_Main_Files' folder)
                for pval = [2, 5, probspecs.n]
                    stodars_option.SubspaceDim = pval;
                    X = stodars_algorithm(funs, X0, stodars_option, probspecs);
                end
                stodars_option.RandomSubspaceMode = 0;
                X = stodars_algorithm(funs, X0, stodars_option, probspecs);
            end
        end
    end
end
