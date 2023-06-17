% Workspace Preparation
% This script is meant to be run at the beginning of each script in this
% project to prepare MATLAB with paths and other code that is redundant in
% each script
%
% Matt Kmiecik
%
% Started 17 JUNE 2023
%

% Sets working directory ----
main_dir = 'M:\luc-stim'; % creating the main dir
cd(main_dir); % setting the working directory

% Directory paths ----
data_dir = fullfile(main_dir, 'data\'); % creating the data folder
output_dir = fullfile(main_dir, 'output\'); % creating the output dir

% Starts EEGLAB ----
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;

% Loads in participant information ----
[NUM,TXT,RAW] = xlsread('doc\ss-info.xlsx');

% Global vars ----
this_elp = fullfile(data_dir, '142130.elp');