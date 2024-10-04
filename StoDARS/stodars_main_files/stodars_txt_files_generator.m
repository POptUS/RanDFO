%% Output files of the StoDARS algorithm
%
%%
% Argonne National Laboratory (USA) / Lawrence Berkeley National Laboratory (USA)
%
% K. Joseph Dzahini (https://github.com/kwassi)
% and Stefan M. Wild (https://wildsm.github.io)
% October 2024
%%

global nprobl fEval_History fEval_Stats ind_sigma

%%
if stodars_option.experiments == 1i % When using stodars_yatsop.m or stodars_benchmark.m
    if stodars_option.FixSeed == 1
        if stodars_option.RandomSubspaceMode == 1
          filestr = ['prob', num2str(nprobl), '_', 'type', probspecs.probtype, '_', ...
            'sigma', num2str(ind_sigma), '_', 'seed', num2str(stodars_option.SeedValue), ...
            '_subdim', num2str(stodars_option.SubspaceDim)];
        else
            filestr = ['prob', num2str(nprobl), '_', 'type', probspecs.probtype, '_', ...
            'sigma', num2str(ind_sigma), '_', 'seed', num2str(stodars_option.SeedValue), ...
            '_fullsp'];
        end
    else
        if stodars_option.RandomSubspaceMode == 1
            filestr = ['prob', num2str(nprobl), '_', 'type', probspecs.probtype, '_', ...
            'sigma', num2str(ind_sigma), '_subdim', num2str(stodars_option.SubspaceDim)];
        else
            filestr = ['prob', num2str(nprobl), '_', 'type', probspecs.probtype, '_', ...
            'sigma', num2str(ind_sigma), '_fullsp'];
        end
    end
    if stodars_option.HistoryFile == 1
        stodars_Hfile = fopen([stodars_option.SaveLocation '/stodars_history_' filestr '.txt'], 'w');
        fprintf(stodars_Hfile, [' %i '  ' %12.20f '  ' %12.20f '  '\n'], fEval_History');
        fclose(stodars_Hfile);
    else
        if stodars_option.StatsFile ~= 1
            stodars_option.SolutionFile = 1;
        end
    end
    if stodars_option.StatsFile == 1
        stodars_AHfile = fopen([stodars_option.SaveLocation '/stodars_stats_' filestr '.txt'], 'w');
        fprintf(stodars_AHfile, [' %i ' repmat(' %12.20f ', 1, probspecs.n) ' %12.20f '  ...
            ' %12.20f '  '\n'], fEval_Stats');
        fclose(stodars_AHfile);
    end
    if stodars_option.SolutionFile == 1
        stodars_Sfile = fopen([stodars_option.SaveLocation '/stodars_solution_' filestr '.txt'], 'w');
        fprintf(stodars_Sfile, [' %12.20f '  '\n'], cur_sol);
        fclose(stodars_Sfile);
    end
else  % When using 'stodars_applications'
    if stodars_option.HistoryFile == 1
        if stodars_option.FixSeed == 1
            filestr = ['history', '_', 'file', '_', 'seed', num2str(stodars_option.SeedValue)];
        else
            filestr = ['history', '_', 'file'];
        end
        stodars_Hfile = fopen([stodars_option.SaveLocation '/stodars_' filestr '.txt'], 'w');
        fprintf(stodars_Hfile, [' %i '  ' %12.20f '  '\n'], fEval_History');
        fclose(stodars_Hfile);
    else
        if stodars_option.StatsFile ~= 1
            stodars_option.SolutionFile = 1;
        end
    end
    if stodars_option.StatsFile == 1
        if stodars_option.FixSeed == 1
            filestr = ['stats', '_', 'file', '_', 'seed', num2str(stodars_option.SeedValue)];
        else
            filestr = ['stats', '_', 'file'];
        end
        stodars_AHfile = fopen([stodars_option.SaveLocation '/stodars_' filestr '.txt'], 'w');
        fprintf(stodars_AHfile, [' %i ' repmat(' %12.20f ', 1, ...
            probspecs.Dimension)  ' %12.20f '  '\n'], fEval_Stats');
        fclose(stodars_AHfile);
    end
    if stodars_option.SolutionFile == 1
        if stodars_option.FixSeed == 1
            filestr = ['solution', '_', 'file', '_', 'seed', num2str(stodars_option.SeedValue)];
        else
            filestr = ['solution', '_', 'file'];
        end
        stodars_Sfile = fopen([stodars_option.SaveLocation '/stodars_' filestr '.txt'], 'w');
        fprintf(stodars_Sfile, [' %12.20f '  '\n'], cur_sol);
        fclose(stodars_Sfile);
    end
end
