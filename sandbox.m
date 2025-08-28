% sandbox.m
% Matt Kmiecik
% place to test functions out

clear; clc;

% Add function path
addpath('fns');

subject_id = '1461831842001_1';
cfg = config();

[success, EEG] = preprocess_eeg(subject_id, cfg);