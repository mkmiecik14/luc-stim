function run_eeg_pipeline(subject_list, options)
% RUN_EEG_PIPELINE - Main orchestration script for modular EEG processing
%
% Usage:
%   run_eeg_pipeline()                    % Process all subjects in config
%   run_eeg_pipeline({'001', '002'})      % Process specific subjects
%   run_eeg_pipeline([], 'resume', true)  % Resume interrupted processing
%
% Inputs:
%   subject_list - Cell array of subject IDs to process (optional)
%   options      - Name-value pairs for processing options:
%                  'resume' - Resume interrupted processing (default: false)
%                  'force_reprocess' - Reprocess existing files (default: false)
%                  'save_intermediate' - Save intermediate results (default: true)
%                  'verbose' - Extra verbose output (default: true)
%
% This script orchestrates the complete EEG processing pipeline:
% 1. Load and preprocess raw EEG data
% 2. Extract resting state blocks  
% 3. Detect artifacts automatically
% 4. Apply ICA (extract-then-ICA approach)
% 5. Interpolate bad channels
% 6. Perform spectral analysis
% 7. Save all results and logs

    % Parse inputs
    if nargin < 1 || isempty(subject_list)
        subject_list = {};  % Will use all subjects from config
    end
    
    if nargin < 2
        options = struct();
    end
    
    % Default options
    if ~isfield(options, 'resume'), options.resume = false; end
    if ~isfield(options, 'force_reprocess'), options.force_reprocess = false; end  
    if ~isfield(options, 'save_intermediate'), options.save_intermediate = true; end
    if ~isfield(options, 'verbose'), options.verbose = true; end
    
    % Initialize pipeline
    fprintf('=== EEG PROCESSING PIPELINE (Modular Version) ===\n');
    fprintf('Started: %s\n\n', datestr(now));
    
    % Add function path
    addpath('fns');
    
    % Load configuration
    try
        fprintf('Loading configuration...\n');
        cfg = config();
    catch ME
        error('Failed to load configuration: %s', ME.message);
    end
    
    % Determine subjects to process
    if isempty(subject_list)
        if isempty(cfg.participants)
            error('No participant information found and no subject list provided');
        end
        subject_list = cfg.participants.ss_list;
        fprintf('Processing all subjects from configuration: %d subjects\n', length(subject_list));
    else
        if ischar(subject_list)
            subject_list = {subject_list};  % Convert single string to cell
        end
        fprintf('Processing specified subjects: %d subjects\n', length(subject_list));
    end
    
    fprintf('Subjects: %s\n\n', strjoin(subject_list, ', '));
    
    % Initialize progress tracking
    n_subjects = length(subject_list);
    processing_start = tic;
    success_count = 0;
    error_count = 0;
    error_log = {};
    
    % Process each subject
    for i = 1:n_subjects
        subject_id = char(subject_list{i});  % Ensure string format
        
        fprintf('\n=== PROCESSING SUBJECT %d/%d: %s ===\n', i, n_subjects, subject_id);
        subject_start = tic;
        
        try
            % Check if already processed (unless force reprocessing)
            if ~options.force_reprocess
                spec_file = fullfile(cfg.paths.output, 'spectral', sprintf('%s-spec-res.mat', subject_id));
                if exist(spec_file, 'file')
                    fprintf('Subject %s already processed. Skipping (use force_reprocess=true to override).\n', subject_id);
                    continue;
                end
            end
            
            % === STEP 1: LOAD RAW EEG ===
            fprintf('\nSTEP 1: Loading raw EEG data...\n');
            EEG_raw = load_raw_eeg(subject_id, cfg);
            
            if options.save_intermediate
                eeg_file = fullfile(cfg.paths.output, 'preprocessed', sprintf('%s-raw.set', subject_id));
                pop_saveset(EEG_raw, 'filename', sprintf('%s-raw.set', subject_id), ...
                           'filepath', fullfile(cfg.paths.output, 'preprocessed'));
            end
            
            % === STEP 2: PREPROCESS EEG ===
            fprintf('\nSTEP 2: Preprocessing EEG data...\n');
            EEG_prepro = preprocess_eeg(EEG_raw, cfg);
            
            if options.save_intermediate
                pop_saveset(EEG_prepro, 'filename', sprintf('%s-prepro.set', subject_id), ...
                           'filepath', fullfile(cfg.paths.output, 'preprocessed'));
            end
            
            % === STEP 3: EXTRACT RESTING STATE BLOCKS ===
            fprintf('\nSTEP 3: Extracting resting state blocks...\n');
            EEG_blocks = extract_blocks(EEG_prepro, cfg);
            
            if options.save_intermediate
                pop_saveset(EEG_blocks, 'filename', sprintf('%s-blocks.set', subject_id), ...
                           'filepath', fullfile(cfg.paths.output, 'preprocessed'));
            end
            
            % === STEP 4: DETECT ARTIFACTS ===
            fprintf('\nSTEP 4: Detecting artifacts...\n');
            [bad_channels, artifact_info] = detect_artifacts(EEG_blocks, cfg);
            
            % === STEP 5: APPLY ICA (with bad channel removal) ===
            fprintf('\nSTEP 5: Applying ICA...\n');
            [EEG_ica_clean, ica_info] = apply_ica(EEG_blocks, bad_channels, cfg);
            
            if options.save_intermediate
                pop_saveset(EEG_ica_clean, 'filename', sprintf('%s-ica-clean.set', subject_id), ...
                           'filepath', fullfile(cfg.paths.output, 'ica'));
            end
            
            % === STEP 6: INTERPOLATE BAD CHANNELS ===
            fprintf('\nSTEP 6: Interpolating bad channels...\n');
            % Need original channel locations from before bad channel removal
            orig_chanlocs = EEG_blocks.chanlocs;  % Channel locations before removal
            EEG_final = interpolate_channels(EEG_ica_clean, bad_channels, orig_chanlocs, cfg);
            
            % === STEP 7: SPECTRAL ANALYSIS ===
            fprintf('\nSTEP 7: Performing spectral analysis...\n');
            spec_results = spectral_analysis(EEG_final, cfg);
            
            % === STEP 8: SAVE ALL RESULTS ===
            fprintf('\nSTEP 8: Saving results...\n');
            save_results(EEG_final, spec_results, artifact_info, ica_info, cfg);
            
            % Success!
            subject_time = toc(subject_start);
            success_count = success_count + 1;
            
            fprintf('\n✓ Subject %s completed successfully in %.1f minutes.\n', ...
                    subject_id, subject_time/60);
            
            % Memory cleanup
            clear EEG_raw EEG_prepro EEG_blocks EEG_ica_clean EEG_final spec_results;
            
        catch ME
            % Handle processing error
            subject_time = toc(subject_start);
            error_count = error_count + 1;
            error_msg = sprintf('Subject %s failed after %.1f minutes: %s', ...
                              subject_id, subject_time/60, ME.message);
            error_log{end+1} = error_msg;
            
            fprintf('\n✗ %s\n', error_msg);
            
            % Save error log
            try
                error_file = fullfile(cfg.paths.output, 'logs', ...
                                    sprintf('%s-error-%s.txt', subject_id, datestr(now, 'yyyymmdd_HHMMSS')));
                fid = fopen(error_file, 'w');
                if fid ~= -1
                    fprintf(fid, 'Error processing subject: %s\n', subject_id);
                    fprintf(fid, 'Timestamp: %s\n', datestr(now));
                    fprintf(fid, 'Error message: %s\n', ME.message);
                    fprintf(fid, 'Stack trace:\n');
                    for j = 1:length(ME.stack)
                        fprintf(fid, '  %s (line %d): %s\n', ME.stack(j).file, ...
                               ME.stack(j).line, ME.stack(j).name);
                    end
                    fclose(fid);
                end
            catch
                % Error logging failed, continue
            end
            
            % Continue with next subject
            continue;
        end
    end
    
    % Final summary
    total_time = toc(processing_start);
    
    fprintf('\n=== PIPELINE PROCESSING COMPLETE ===\n');
    fprintf('Total time: %.1f minutes (%.1f hours)\n', total_time/60, total_time/3600);
    fprintf('Subjects processed: %d\n', n_subjects);
    fprintf('Successful: %d (%.1f%%)\n', success_count, 100*success_count/n_subjects);
    fprintf('Failed: %d (%.1f%%)\n', error_count, 100*error_count/n_subjects);
    
    if success_count > 0
        fprintf('Average time per subject: %.1f minutes\n', (total_time/60)/success_count);
    end
    
    % Report errors
    if error_count > 0
        fprintf('\nERRORS ENCOUNTERED:\n');
        for i = 1:length(error_log)
            fprintf('  %d. %s\n', i, error_log{i});
        end
        fprintf('Check log files in %s for detailed error information.\n', ...
                fullfile(cfg.paths.output, 'logs'));
    else
        fprintf('\n✓ All subjects processed successfully!\n');
    end
    
    % Final output summary
    fprintf('\nOUTPUT SUMMARY:\n');
    fprintf('  Preprocessed EEG: %s/preprocessed/\n', cfg.paths.output);
    fprintf('  ICA results: %s/ica/\n', cfg.paths.output);
    fprintf('  Spectral results: %s/spectral/\n', cfg.paths.output);
    fprintf('  Processing logs: %s/logs/\n', cfg.paths.output);
    fprintf('  Quality control: %s/qc/\n', cfg.paths.output);
    
    fprintf('\nResults are compatible with original pipeline format.\n');
    fprintf('Pipeline completed: %s\n', datestr(now));
    
end