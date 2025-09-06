% BATCH_SPECTRAL - Batch spectral analysis of EEG data for all subjects
% 
% This script processes all subjects listed in cfg.participants.ss_list
% using the spectral_analysis() function. Reports any failed subjects at completion.

clear; clc;

fprintf('Starting batch spectral analysis...\n\n');

% Load configuration
cfg = config();

% Get subject list
subjects = cfg.participants.ss_list;
fprintf('Found %d subjects to process\n\n', length(subjects));

% Initialize tracking
failed_subjects = {};

% Process each subject
for i = 1:length(subjects)
    subject_id = subjects{i};
    fprintf('Processing subject %d/%d: %s\n', i, length(subjects), subject_id);
    
    [success, ~] = spectral_analysis(subject_id, cfg);
    
    if ~success
        failed_subjects{end+1} = subject_id;
        fprintf('  FAILED: Subject %s\n', subject_id);
    else
        fprintf('  SUCCESS: Subject %s\n', subject_id);
    end
    fprintf('\n');
end

% Final summary
fprintf('=== BATCH SPECTRAL ANALYSIS SUMMARY ===\n');
fprintf('Total subjects: %d\n', length(subjects));
fprintf('Successful: %d\n', length(subjects) - length(failed_subjects));
fprintf('Failed: %d\n', length(failed_subjects));

if ~isempty(failed_subjects)
    fprintf('\nFailed subjects: %s\n', strjoin(failed_subjects, ', '));
else
    fprintf('\nAll subjects processed successfully!\n');
end