function EEG = load_raw_eeg(subject_id, cfg)
% LOAD_RAW_EEG - Load and setup raw EEG data from BDF file
%
% Usage:
%   EEG = load_raw_eeg(subject_id, cfg)
%
% Inputs:
%   subject_id - String or char, subject identifier (e.g., '001')
%   cfg        - Configuration structure from config()
%
% Outputs:
%   EEG        - EEGLAB EEG structure with loaded and configured data
%
% This function:
% 1. Loads BDF file using pop_biosig
% 2. Removes external channels
% 3. Keeps only 64-channel Biosemi montage
% 4. Configures channel locations
% 5. Sets reference configuration

    fprintf('Loading raw EEG data for subject %s...\n', subject_id);
    
    % Find BDF file first
    bdf_pattern = fullfile(cfg.paths.data, sprintf('%s.bdf', subject_id));
    bdf_files = dir(bdf_pattern);

    if ~isempty(bdf_files)
        % BDF file found - use existing logic
        if length(bdf_files) > 1
            error('Multiple BDF files found for subject %s. Expected exactly one file matching pattern: %s', subject_id, bdf_pattern);
        end

        bdf_file = fullfile(bdf_files(1).folder, bdf_files(1).name);
        fprintf('Found BDF file: %s\n', bdf_files(1).name);

        % Load raw data using biosemi import
        try
            EEG = pop_biosig(...
                bdf_file,...
                'ref', [1], ...
                'refoptions', {'keepref','on'},...
                'importannot', 'off',... % does not import EDF annotations
                'bdfeventmode', 6 ...     % event mode that syncs with EMSE events
            );
            fprintf('Successfully loaded BDF file: %d channels, %d samples, %.1f Hz\n', ...
                    EEG.nbchan, EEG.pnts, EEG.srate);
        catch ME
            error('Failed to load BDF file: %s', ME.message);
        end

    else
        % BDF not found - check for stitched file in output/raw-data
        stitched_pattern = fullfile(cfg.paths.output, 'raw-data', sprintf('%s.set', subject_id));
        stitched_files = dir(stitched_pattern);

        if ~isempty(stitched_files)
            % Stitched file found
            if length(stitched_files) > 1
                error('Multiple stitched files found for subject %s. Expected exactly one file matching pattern: %s', subject_id, stitched_pattern);
            end

            stitched_file = fullfile(stitched_files(1).folder, stitched_files(1).name);
            fprintf('BDF not found, using stitched file: %s\n', stitched_files(1).name);

            % Load stitched data using EEGLAB import
            try
                EEG = pop_loadset('filename', stitched_files(1).name, ...
                                 'filepath', stitched_files(1).folder);
                fprintf('Successfully loaded stitched file: %d channels, %d samples, %.1f Hz\n', ...
                        EEG.nbchan, EEG.pnts, EEG.srate);
            catch ME
                error('Failed to load stitched file: %s', ME.message);
            end

        else
            % Neither BDF nor stitched file found
            error('No BDF file or stitched file found for subject %s. Check %s and %s', ...
                  subject_id, cfg.paths.data, fullfile(cfg.paths.output, 'raw-data'));
        end
    end
    
    % Remove external channels that are not being used
    if any(ismember(cfg.channels.external_remove, {EEG.chanlocs.labels}))
        fprintf('Removing external channels...\n');
        EEG = pop_select(EEG, 'rmchannel', cfg.channels.external_remove);
    end
    
    % Check if more than 64 channels and keep only A1-A32, B1-B32
    if EEG.nbchan > 64
        fprintf('More than 64 channels detected (%d). Keeping only A1-A32, B1-B32...\n', EEG.nbchan);
        % Find which channels from our keep list actually exist
        existing_chans = {EEG.chanlocs.labels};
        chans_to_keep = cfg.channels.keep(ismember(cfg.channels.keep, existing_chans));
        EEG = pop_select(EEG, 'channel', chans_to_keep);
    else
        fprintf('64 or fewer channels detected (%d); not removing any channels.\n', EEG.nbchan);
    end
    
    % Configure channel locations
    if exist(cfg.channels.chan_info, 'file') && exist(cfg.channels.chan_locs, 'file')
        fprintf('Loading channel locations from: %s\n', cfg.paths.doc);
        try
            fprintf('  Configuring channel locations...\n');
            load(cfg.channels.chan_info); % channel info
            load(cfg.channels.chan_locs); % channel locs
            EEG.chaninfo = chan_info;
            EEG.chanlocs = chan_locs;
            
            fprintf('Successfully configured channel locations for %d channels.\n', EEG.nbchan);
        catch ME
            warning(sprintf('Could not load channel locations: %s', ME.message));
        end
    else
        warning('Channel location files not found');
    end
    
    % Set reference configuration (A1/Fp1 as reference, as per original)
    try
        EEG = pop_chanedit(EEG, 'setref', {'1:64' 'Fp1'}); 
        fprintf('Set reference configuration.\n');
    catch ME
        warning(sprintf('Could not set reference: %s', ME.message));
    end
    
    % Add processing history
    EEG.etc.processing_history = {};
    EEG.etc.processing_history{1} = sprintf('Loaded raw data: %s', datetime('now'));
    EEG.etc.processing_history{2} = sprintf('Original channels: %d, Final channels: %d', ...
                                          length(bdf_files), EEG.nbchan);
    
    % Update dataset name
    EEG.setname = sprintf('%s-raw', subject_id);
    EEG = eeg_checkset(EEG);
    
    fprintf('Raw EEG loading completed for subject %s.\n', subject_id);
    
end