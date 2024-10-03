% This script provides to stars_yatsop.m, stars_applications.m and
% stars_benchmark the paths for the main algorithm and problems files,
% which are necessary for each run of the STARS algorithm.
% It also creates a folder for the outputs txt files.
%
% K. Joseph Dzahini (https://github.com/kwassi)
% and Stefan M. Wild (https://wildsm.github.io)
%%
% This code and its updates are available at https://github.com/POptUS/RanDFO
%%
%  Argonne National Laboratory, Mathematics and Computer Science Division
%      Kwassi Joseph Dzahini and Stefan M. Wild, October 2022
%%

% Main files location
main_files_location = 'stars_main_files';

% Adding main files location folder to the MATLAB path.
addpath([main_files_location, '/']);
addpath([main_files_location, '/CSV_DFO/']);
addpath([main_files_location, '/Pounders_files/']);
addpath([main_files_location, '/model_constructor/']);

% Locations of the two dependencies (loaded as submodules)
bendfo_location = 'problems/BenDFO';
yatsop_location = 'problems/YATSOp';

% Adding YATSOp folder (containing experiment main files) to the MATLAB path.
if stars_option.yatsop == 1
    addpath([yatsop_location, '/m/']);
end

% Adding benchmark folder to the MATLAB path.
if stars_option.benchmark == 1
    addpath([bendfo_location, '/m/']);
    addpath([bendfo_location, '/data/']);
    addpath('problems/benchmark/m/'); % For problems_specifications.m
end

stars_default_options;

%% Location for the outputs (see txt_files_generator.m)
stars_option.SaveLocation = 'stars_outputs';
if ~exist('stars_outputs', 'dir')
    mkdir stars_outputs;
end
