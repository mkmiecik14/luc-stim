
function [success, EEG] = apply_ica(subject_id, cfg)
% APPLY_ICA - Remove bad channels, apply ICA, and reject artifactual components
%
% Usage:
%   [success, EEG] = apply_ica(subject_id, cfg)
%
% Inputs:
%   subject_id - String or char, subject identifier (e.g., '001')
%   cfg        - Configuration structure from config()
%
% Outputs:
%   success - Logical flag indicating successful ICA processing
%   EEG     - EEGLAB EEG structure with bad channels removed, ICA applied, and ICs rejected (empty if success = false)
%
% This function implements the "remove → ICA → clean" approach:
% 1. Load preprocessed data
% 2. Detect and remove bad channels before ICA
% 3. Run ICA decomposition on clean channel set
% 4. Classify and reject artifactual independent components
% 5. Save results to output/ica directory

    fprintf('Applying ICA to subject %s...\n', subject_id);
    
    % Initialize outputs
    success = false;
    EEG = [];
    
    % Load preprocessed data
    prepro_file = fullfile(cfg.paths.output, 'preprocessed', sprintf('%s-prepro.set', subject_id));
    if ~exist(prepro_file, 'file')
        fprintf('ERROR: Preprocessed file not found for subject %s: %s\n', subject_id, prepro_file);
        return;
    end
    
    try
        EEG = pop_loadset('filename', sprintf('%s-prepro.set', subject_id), ...
                          'filepath', fullfile(cfg.paths.output, 'preprocessed'));
        fprintf('Successfully loaded preprocessed data for subject %s\n', subject_id);
    catch ME
        fprintf('ERROR: Failed to load preprocessed data for subject %s: %s\n', subject_id, ME.message);
        return;
    end
    
    % Store original parameters
    orig_n_chans = EEG.nbchan;
    
    % Begin ICA processing with comprehensive error handling
    try
        % Detect and remove bad channels
        fprintf('Detecting bad channels...\n');
        orig_labels = {EEG.chanlocs.labels};
        
        EEG = pop_clean_rawdata(EEG, ...
            'FlatlineCriterion', cfg.params.bad_chan.FlatlineCriterion, ...     % Channel: flat channels
            'ChannelCriterion', cfg.params.bad_chan.ChannelCriterion, ...       % Channel: correlation-based
            'LineNoiseCriterion', cfg.params.bad_chan.LineNoiseCriterion, ...   % Channel: line noise
            'Highpass', 'off', ...                % No filtering
            'BurstCriterion', 'off', ...          % No burst detection  
            'WindowCriterion', 'off', ...         % No window rejection
            'BurstRejection', 'off', ...          % No time rejection
            'Distance', 'Euclidian');
        
        % Determine which channels were removed
        remaining_labels = {EEG.chanlocs.labels};
        removed_channels = setdiff(orig_labels, remaining_labels);
        
        if ~isempty(removed_channels)
            fprintf('Removed %d bad channels: %s\n', length(removed_channels), strjoin(removed_channels, ', '));
        else
            fprintf('No bad channels detected.\n');
        end
    
        fprintf('Channels after bad channel removal: %d/%d (%.1f%% retained)\n', ...
                EEG.nbchan, orig_n_chans, 100 * EEG.nbchan / orig_n_chans);
    
        % Calculate rank for ICA
        % After linked mastoid referencing and bad channel removal, 
        % the rank is typically n_channels - 1
        ica_rank = orig_n_chans - 1 - length(removed_channels);
        
        fprintf('Using ICA rank: %d (channels: %d)\n', ica_rank, EEG.nbchan);
    
        % Run ICA decomposition
        fprintf('Running ICA decomposition...\n');
        EEG = pop_runamica(EEG, 'pcakeep', ica_rank);
        
        %EEG = pop_runica(EEG,'icatype', 'runica','extended',1, ...
        %    'interrupt','on','pca',ica_rank);
        
        fprintf('ICA decomposition completed: %d components\n', size(EEG.icaweights, 1));
        
        % 7. Update dataset name and save results
        EEG.setname = sprintf('%s-ica', subject_id);
        EEG = eeg_checkset(EEG);
        
        % Save ICA-cleaned data to output/ica directory
        fprintf('Saving ICA-cleaned data for subject %s...\n', subject_id);
        filename = sprintf('%s-ica.set', subject_id);
        EEG = pop_saveset(EEG, 'filename', filename, 'filepath', fullfile(cfg.paths.output, 'ica'));
        
        % Set success flag
        success = true;
        fprintf('ICA processing completed for subject %s:\n', subject_id);
        fprintf('  Original channels: %d\n', orig_n_chans);
        fprintf('  Final channels: %d\n', EEG.nbchan);
        fprintf('  ICA components: %d (rank: %d)\n', size(EEG.icaweights, 1), ica_rank);
        
    catch ME
        fprintf('ERROR: ICA processing failed for subject %s: %s\n', subject_id, ME.message);
        if ~isempty(ME.stack)
            fprintf('Error occurred at: %s\n', ME.stack(1).name);
        end
        success = false;
        EEG = [];
        return;
    end

end