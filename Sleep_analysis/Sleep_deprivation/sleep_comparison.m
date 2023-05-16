%% 1. sleep comparison
% read in the control csv files
% binned activity control
tic
[filename_csv, PathName, FilterIndex]=uigetfile('*');
cd(PathName);
fopen(filename_csv);             
binned_act_ctrl = readtable(filename_csv,'PreserveVariableNames',true)
toc
%% read in the csv files
% binned rest control
tic
[filename_csv, PathName, FilterIndex]=uigetfile('*');
cd(PathName);
fopen(filename_csv);             
binned_rest_ctrl = readtable(filename_csv,'PreserveVariableNames',true)
toc
%% read in the csv files
% binned Wake activity control
tic
[filename_csv, PathName, FilterIndex]=uigetfile('*');
cd(PathName);
fopen(filename_csv);             
binned_wakeAct_ctrl = readtable(filename_csv,'PreserveVariableNames',true)
toc
%% read in the csv files
% cumulative rest day 1 control
tic
[filename_csv, PathName, FilterIndex]=uigetfile('*');
cd(PathName);
fopen(filename_csv);             
cRest_day1_ctrl = readtable(filename_csv,'PreserveVariableNames',true)
toc
%% read in the csv files
% cumulative Rest Night1 control
tic
[filename_csv, PathName, FilterIndex]=uigetfile('*');
cd(PathName);
fopen(filename_csv);             
cRest_night1_ctrl = readtable(filename_csv,'PreserveVariableNames',true)
toc
%% read in the csv files
% cumulative rest testphase
tic
[filename_csv, PathName, FilterIndex]=uigetfile('*');
cd(PathName);
fopen(filename_csv);             
cRest_testphase_ctrl = readtable(filename_csv,'PreserveVariableNames',true)
toc
%% read in the csv files
% rest bout length
tic
[filename_csv, PathName, FilterIndex]=uigetfile('*');
cd(PathName);
fopen(filename_csv);             
restBout_length_ctrl = readtable(filename_csv,'PreserveVariableNames',true)
toc
%% read in the csv files
% rest bout number control
tic
[filename_csv, PathName, FilterIndex]=uigetfile('*');
cd(PathName);
fopen(filename_csv);             
restBout_number_ctrl = readtable(filename_csv,'PreserveVariableNames',true)
toc
%% 2. treated fish (could be either sleep-deprived or tumor-induced,...)
% read in the csv files from the treated group of the **same** day
% binned activity treatment
tic
[filename_csv, PathName, FilterIndex]=uigetfile('*');
cd(PathName);
fopen(filename_csv);             
binned_act_treatment = readtable(filename_csv,'PreserveVariableNames',true)
toc
%% binned rest treatment
tic
[filename_csv, PathName, FilterIndex]=uigetfile('*');
cd(PathName);
fopen(filename_csv);             
binned_rest_treatment = readtable(filename_csv,'PreserveVariableNames',true)
toc
%% binned wake activity treatment
tic
[filename_csv, PathName, FilterIndex]=uigetfile('*');
cd(PathName);
fopen(filename_csv);             
binned_wakeAct_treatment = readtable(filename_csv,'PreserveVariableNames',true)
toc
%% cumulative rest day1 treatment
tic
[filename_csv, PathName, FilterIndex]=uigetfile('*');
cd(PathName);
fopen(filename_csv);             
cRest_treatment = readtable(filename_csv,'PreserveVariableNames',true)
toc
%% cumulative rest night1 treatment
tic
[filename_csv, PathName, FilterIndex]=uigetfile('*');
cd(PathName);
fopen(filename_csv);             
cRest_night1_treatment = readtable(filename_csv,'PreserveVariableNames',true)
toc
%% cumulative rest testphase treatment
tic
[filename_csv, PathName, FilterIndex]=uigetfile('*');
cd(PathName);
fopen(filename_csv);             
cRest_test_treatment = readtable(filename_csv,'PreserveVariableNames',true)
toc
%% restbout length treatment
tic
[filename_csv, PathName, FilterIndex]=uigetfile('*');
cd(PathName);
fopen(filename_csv);             
restBout_length_treatment = readtable(filename_csv,'PreserveVariableNames',true)
toc
%% restbout number treatment
tic
[filename_csv, PathName, FilterIndex]=uigetfile('*');
cd(PathName);
fopen(filename_csv);             
restBout_number_treatment = readtable(filename_csv,'PreserveVariableNames',true)
toc
%% 3.Plots
%% 9. ACTIVITY/REST PLOTS
% this calculates:
% waking activity (seconds of motion per 'waking minute')
%   where a 'waking minute' is a minute with any activity
% rest (number of resting minutes in each 10 minute interval)
% activity (number of active seconds per 10 minute interval)

geno = {'ctrl'; ctrl};

%To plot the designated timeframe based on d96_clip:
% Save as svg for best results in Illustrator, etc
% MANUALLY SAVE from window to preserve display dimensions
ctrl = binned_wakeAct_ctrl
treatment = binned_wakeAct_treatment
f = fct_crSleepPlots_comparison(ctrl,treatment,geno,expName);

% Export structure f by genotype in case you want to replot outside Matlab
T_f = struct2table(f, 'AsArray', 1);
csvwrite('output_binned_rest_ctrl.csv', T_f.rest{1}); % ctrl
csvwrite('output_binned_wakeAct_ctrl.csv', T_f.wakeAct{1}); % ctrl
csvwrite('output_binned_act_ctrl.csv', T_f.activity{1}); % ctrl



