% This function computes function estimates using the average of the noisy
% observations. It also indicates to the main STARS algorithm whether the
% maximum function evaluation is reached or not.
% When the maximum evaluations is reached, the function returns the complex
% number '1i' (i.e., the complex imaginary number i) which is not stored in
% the history files but simply indicates to the algorithm to stop.
%
% This function, together with 'txt_files_generator.m', also handles the
% generation of .txt Solutions, History and Stats files.
%
%% Below:
%
% nfEval = Number of function evaluations
%
% When using stars_yatsop.m or stars_benchmark.m,
% fEval_History = a matrix with rows [#Funct_eval, smooth_funct_value,
% noisy_value] (for the history file);
% fEval_Stats = [#Funct_eval, x, smooth_funct_value, noisy_value] (for the stats file)
%
% When using stars_applications.m
% fEval_History = a matrix with rows [#Funct_eval, noisy_value]
% fEval_Stats = [#Funct_eval, x, noisy_value]
%
% K. Joseph Dzahini (https://github.com/kwassi)
% and Stefan M. Wild (https://wildsm.github.io)
%%
% This code and its updates are available at https://github.com/POptUS/RanDFO
%%
%  Argonne National Laboratory, Mathematics and Computer Science Division
%      Kwassi Joseph Dzahini and Stefan M. Wild, October 2022
%%

function z = MCestimate(x, sample_size, funs, stars_option)

my_fun = funs.myfun;

global nfEval fEval_History fEval_Stats nfEvalExceeded

y = zeros(sample_size, 1);
if stars_option.experiments == 1i
    smooth_fun = funs.mysmoothfun;
end
for i = 1:sample_size
    if nfEval <= stars_option.MaxFuncEval
        obj = my_fun(x);
        if length(obj) > 1
            warning('Error defining the objective function.');
            warning('Outputs of the objective function should be scalars');
            nfEvalExceeded = 1;
            break
        else
            y(i) = obj;
        end
        if stars_option.HistoryFile == 1
            if stars_option.experiments == 1i
                fEval_History(nfEval, :) = [nfEval, smooth_fun(x), y(i)];
            else
                fEval_History(nfEval, :) = [nfEval, y(i)];  % Before any update, see
                %                                             txt_files_generator.m
            end
        end
        if stars_option.StatsFile == 1
            if iscolumn(x)
                x = x';
            end
            if stars_option.experiments == 1i
                fEval_Stats(nfEval, :) = [nfEval, x, smooth_fun(x), y(i)];
            else
                fEval_Stats(nfEval, :) = [nfEval, x, y(i)];
            end
        end
        nfEval = nfEval + 1;
    else
        nfEvalExceeded = 1;
        break
    end
end
if nfEvalExceeded ~= 1
    z = mean(y);
else
    z = 1i;
end
end
