function [success, EEG] = preprocess_eeg(subject_id, cfg)
% PREPROCESS_EEG - Apply basic preprocessing to raw EEG data
%
% Usage:
%   [success, EEG] = preprocess_eeg(subject_id, cfg)
%
% Inputs:
%   subject_id - String or char, subject identifier (e.g., '001')
%   cfg - Configuration structure from config()
%
% Outputs:
%   success - Logical flag indicating successful preprocessing
%   EEG - Preprocessed EEGLAB EEG structure (empty if success = false)
%
% This function applies:
% 1. Downsampling to target sampling rate
% 2. Common average reference
% 3. DC offset removal
% 4. High-pass filtering
% 5. Line noise removal (cleanline)

    fprintf('Preprocessing EEG data for subject %s...\n', subject_id);
    
    % Initialize success flag
    success = false;
    EEG = [];
    
    % Load raw EEG data
    try
        EEG = load_raw_eeg(subject_id, cfg);
        fprintf('Successfully loaded raw EEG data for subject %s\n', subject_id);
    catch ME
        fprintf('ERROR: Failed to load raw EEG data for subject %s: %s\n', subject_id, ME.message);
        return;
    end
    
    % Store original parameters
    orig_srate = EEG.srate;
    orig_pnts = EEG.pnts;
    
    % Begin preprocessing steps with comprehensive error handling
    try
        % 1. Downsample to target sampling rate
        if EEG.srate ~= cfg.params.srate_target
            fprintf('Downsampling from %.1f Hz to %d Hz...\n', EEG.srate, cfg.params.srate_target);
            EEG = pop_resample(EEG, cfg.params.srate_target);
            fprintf('Downsampling completed: %d samples -> %d samples\n', orig_pnts, EEG.pnts);
        else
            fprintf('Sampling rate already at target (%.1f Hz), skipping downsampling.\n', EEG.srate);
        end
    
        % 2. Re-reference to linked mastoids
        fprintf('  Re-referencing to linked mastoids...\n');
        EEG = pop_reref(EEG, cfg.params.linked_mast_ref, 'keepref', 'on');
    
        % 4. High-pass filter
        fprintf('Applying high-pass filter at %.1f Hz...\n', cfg.params.highpass_freq);
        % Using the same parameters as original: -6dB @ 2Hz, transition band width
        EEG = pop_eegfiltnew(EEG, 'locutoff', cfg.params.highpass_freq, 'plotfreqz', 0);
        fprintf('High-pass filtering completed.\n');
    
        % 5. Remove electrical line noise using cleanline
        fprintf('Removing %d Hz line noise using cleanline...\n', cfg.params.line_freq);
        EEG = pop_cleanline(EEG, ...
            'bandwidth', 2, ...
            'chanlist', [1:EEG.nbchan], ...
            'computepower', 1, ...
            'linefreqs', cfg.params.line_freq, ...
            'normSpectrum', 0, ...
            'p', 0.01, ...
            'pad', 2, ...
            'plotfigures', 0, ...
            'scanforlines', 1, ...
            'sigtype', 'Channels', ...
            'tau', 100, ...
            'verb', 1, ...
            'winsize', 4, ...
            'winstep', 1);
        
        % Extracting blocks
        fprintf('Extracting resting state blocks...');
        EEG = extract_blocks(EEG, cfg);

        % Update dataset name
        EEG.setname = sprintf('%s-prepro', subject_id);
        EEG = eeg_checkset(EEG);

        % Saving out participant preprocessed data
        fprintf('Saving out preprocessed data for subject %s:\n', subject_id);
        filename = [EEG.setname '.set'];
        EEG = pop_saveset(EEG, 'filename', filename, 'filepath', fullfile(cfg.paths.output, 'preprocessed'));
        
        % Set success flag if we reach this point
        fprintf('Preprocessing completed for subject %s:\n', subject_id);
        success = true;
        
    catch ME
        fprintf('ERROR: Preprocessing failed for subject %s: %s\n', subject_id, ME.message);
        fprintf('Error occurred at: %s\n', ME.stack(1).name);
        success = false;
        EEG = [];
        return;
    end
    
end