%% Output files of the STARS algorithm
%
% K. Joseph Dzahini (https://github.com/kwassi)
% and Stefan M. Wild (https://wildsm.github.io)
%%
% This code and its updates are available at https://github.com/POptUS/RanDFO
%%
%  Argonne National Laboratory, Mathematics and Computer Science Division
%      Kwassi Joseph Dzahini and Stefan M. Wild, October 2022
%%

global nprobl fEval_History fEval_Stats ind_sigma

%%
if stars_option.experiments == 1i % When using stars_yatsop.m or stars_benchmark.m
    if stars_option.FixSeed == 1
        filestr = ['prob', num2str(nprobl), '_', 'type', probspecs.probtype, '_', 'sigma', ...
            num2str(ind_sigma), '_', 'p', num2str(stars_option.SubspaceDim), '_', ...
            'seed', num2str(stars_option.SeedValue)];
    else
        filestr = ['prob', num2str(nprobl), '_', 'type', probspecs.probtype, '_', 'sigma', ...
            num2str(ind_sigma), '_', 'p', num2str(stars_option.SubspaceDim)];
    end
    if stars_option.HistoryFile == 1
        stars_Hfile = fopen([stars_option.SaveLocation '/history_' filestr '.txt'], 'w');
        fprintf(stars_Hfile, [' %i '  ' %12.20f '  ' %12.20f '  '\n'], fEval_History');
        fclose(stars_Hfile);
    else
        if stars_option.StatsFile ~= 1
            stars_option.SolutionFile = 1;
        end
    end
    if stars_option.StatsFile == 1
        stars_AHfile = fopen([stars_option.SaveLocation '/stats_' filestr '.txt'], 'w');
        fprintf(stars_AHfile, [' %i ' repmat(' %12.20f ', 1, probspecs.n) ' %12.20f '  ...
            ' %12.20f '  '\n'], fEval_Stats');
        fclose(stars_AHfile);
    end
    if stars_option.SolutionFile == 1
        stars_Sfile = fopen([stars_option.SaveLocation '/solution_' filestr '.txt'], 'w');
        fprintf(stars_Sfile, [' %12.20f '  '\n'], cur_sol);
        fclose(stars_Sfile);
    end
else  % When using 'stars_applications'
    if stars_option.HistoryFile == 1
        if stars_option.FixSeed == 1
            filestr = ['history', '_', 'file', '_', 'p', num2str(stars_option.SubspaceDim), '_', ...
                'seed', num2str(stars_option.SeedValue)];
        else
            filestr = ['history', '_', 'file', '_', 'p', num2str(stars_option.SubspaceDim)];
        end
        stars_Hfile = fopen([stars_option.SaveLocation '/stars_' filestr '.txt'], 'w');
        fprintf(stars_Hfile, [' %i '  ' %12.20f '  '\n'], fEval_History');
        fclose(stars_Hfile);
    else
        if stars_option.StatsFile ~= 1
            stars_option.SolutionFile = 1;
        end
    end
    if stars_option.StatsFile == 1
        if stars_option.FixSeed == 1
            filestr = ['stats', '_', 'file', '_', 'p', num2str(stars_option.SubspaceDim), '_', ...
                'seed', num2str(stars_option.SeedValue)];
        else
            filestr = ['stats', '_', 'file', '_', 'p', num2str(stars_option.SubspaceDim)];
        end
        stars_AHfile = fopen([stars_option.SaveLocation '/stars_' filestr '.txt'], 'w');
        fprintf(stars_AHfile, [' %i ' repmat(' %12.20f ', 1, ...
            probspecs.Dimension)  ' %12.20f '  '\n'], fEval_Stats');
        fclose(stars_AHfile);
    end
    if stars_option.SolutionFile == 1
        if stars_option.FixSeed == 1
            filestr = ['file', '_', 'p', num2str(stars_option.SubspaceDim), ...
                '_', 'seed', num2str(stars_option.SeedValue)];
        else
            filestr = ['file', '_', 'p', num2str(stars_option.SubspaceDim)];
        end
        stars_Sfile = fopen([stars_option.SaveLocation '/stars_solution_' filestr '.txt'], 'w');
        fprintf(stars_Sfile, [' %12.20f '  '\n'], cur_sol);
        fclose(stars_Sfile);
    end
end
