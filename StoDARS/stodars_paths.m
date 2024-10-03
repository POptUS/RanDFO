% This script provides to stodars_yatsop.m, stodars_benchmark.m and
% stodars_applications.m, the paths for the main algorithm and problems
% files, which are necessary for each run of the stodars algorithm.
% It also creates a folder for the outputs txt files.
%
%%
% Argonne National Laboratory (USA) / Lawrence Berkeley National Laboratory (USA)
%
% K. Joseph Dzahini (https://github.com/kwassi)
% and Stefan M. Wild (https://wildsm.github.io)
% October 2024
%
%%
% Main files location
main_files_location = 'stodars_main_files';

% Adding main files location folder to the MATLAB path.
addpath([main_files_location, '/']);

% Locations of the two dependencies (loaded as submodules)
bendfo_location = 'problems/BenDFO';
yatsop_location = 'problems/YATSOp';

% Adding YATSOp folder (containing experiment main files) to the MATLAB path.
if stodars_option.yatsop == 1
    addpath([yatsop_location, '/m/']);
end

% Adding benchmark folder to the MATLAB path.
if stodars_option.benchmark == 1
    addpath([bendfo_location, '/m/']);
    addpath([bendfo_location, '/data/']);
    addpath('problems/benchmark/m/'); % For problems_specifications.m
end

stodars_default_options;

%% Location for the outputs (see txt_files_generator.m)
stodars_option.SaveLocation = 'stodars_outputs';
if ~exist('stodars_outputs', 'dir')
    mkdir stodars_outputs;
end
