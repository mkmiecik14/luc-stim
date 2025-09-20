% run_stitch.m
% Script to stitch BDF files for the problematic subject

% Load configuration
cfg = config();

% Define the subject and their two BDF files
subject_id = '1461831842050';       % subject ID
part1_file = '1461832050_3.bdf';    % before stim
part2_file = '1461842050_3.bdf';    % after stim
outname = '1461831842050_3';        % output name 

% Run the stitching
try
    stitch_bdfs(subject_id, part1_file, part2_file, outname, cfg);
    fprintf('\nStitching completed successfully!\n');
    fprintf('You can now run the normal pipeline on subject %s\n', subject_id);
catch ME
    fprintf('ERROR: Stitching failed: %s\n', ME.message);
end