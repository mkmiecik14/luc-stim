function [success, spec_results] = spectral_analysis(subject_id, cfg)
% SPECTRAL_ANALYSIS - Perform spectral decomposition and IAF analysis
%
% Usage:
%   [success, spec_results] = spectral_analysis(subject_id, cfg)
%
% Inputs:
%   subject_id - String or char, subject identifier (e.g., '001')
%   cfg        - Configuration structure from config()
%
% Outputs:
%   success      - Logical flag indicating successful spectral analysis
%   spec_results - Structure containing spectral analysis results (empty if success = false)
%
% This function performs spectral analysis on individual resting state blocks:
% 1. Load ICA data from output/ica directory
% 2. Classify and reject artifactual independent components
% 3. Spherical interpolation to restore 64-channel montage
% 4. Re-reference to average signal
% 5. Epoch data by blocks (eyes open/closed conditions)
% 6. Calculate power spectral density (PSD) using spectopo
% 7. Compute Individual Alpha Frequency (IAF) using restingIAF
% 8. Save results to output/spectral directory
%
% Note: Pink/white noise correction is NOT applied in this simplified version

    fprintf('Performing spectral analysis for subject %s...\n', subject_id);
    
    % Initialize outputs
    success = false;
    spec_results = [];
    
    % Load ICA data
    ica_file = fullfile(cfg.paths.output, 'ica', sprintf('%s-ica.set', subject_id));
    if ~exist(ica_file, 'file')
        fprintf('ERROR: ICA file not found for subject %s: %s\n', subject_id, ica_file);
        return;
    end
    
    try
        EEG = pop_loadset('filename', sprintf('%s-ica.set', subject_id), ...
                          'filepath', fullfile(cfg.paths.output, 'ica'));
        fprintf('Successfully loaded ICA data for subject %s\n', subject_id);
    catch ME
        fprintf('ERROR: Failed to load ICA data for subject %s: %s\n', subject_id, ME.message);
        return;
    end
    
    % Get spectral analysis parameters
    wsize = cfg.params.wsize;              % Window size in seconds
    overlap = cfg.params.overlap;          % Overlap in seconds  
    frange = cfg.params.spectral.frange;   % Frequency range [Hz]
    alpha_window = cfg.params.spectral.alpha_window; % Alpha search window [Hz]
    cmin = cfg.params.spectral.cmin;       % Min channels for cross-channel avg
    fw = cfg.params.spectral.fw;           % SGF frame width
    poly = cfg.params.spectral.poly;       % SGF polynomial order
    
    % Block definitions
    block_names = cfg.blocks.names;
    block_descriptions = cfg.blocks.descriptions;
    
    fprintf('Spectral parameters: window=%.1fs, overlap=%.1fs, frange=[%.1f-%.1f]Hz, alpha=[%.1f-%.1f]Hz\n', ...
            wsize, overlap, frange(1), frange(2), alpha_window(1), alpha_window(2));
    
    % Initialize output structure
    spec_results = struct();
    spec_results.subject_id = subject_id;
    %spec_results.timestamp = datetime('now');
    spec_results.param_wsize = wsize;
    spec_results.param_overlap = overlap;
    spec_results.param_frange = frange;
    spec_results.param_alpha_window = alpha_window;
    spec_results.param_cmin = cmin;
    spec_results.param_fw = fw;
    spec_results.param_poly = poly;
    spec_results.block_names = block_names;
    spec_results.block_descriptions = block_descriptions;
    
    % Store original channel count before processing
    orig_n_chans = EEG.nbchan;
    
    % Pre-allocate result matrices
    n_blocks = length(block_names);
    
    % Calculate expected frequency resolution
    N = EEG.srate * wsize;          % Number of data points in FFT window
    n_freqs = floor(N/2) + 1;       % Number of frequency bins
    
    % Initialize matrices (will be updated after channel interpolation)
    iaf_summary = NaN(n_blocks, 2); % Summary IAF [PAF, COG] per block

    % Begin spectral analysis processing with comprehensive error handling
    try
        % 1. Classify independent components
        fprintf('Classifying independent components using ICLabel...\n');
        
        try
            EEG = pop_iclabel(EEG, 'default');
            
            % ICLabel classifications: [Brain Muscle Eye Heart LineNoise ChannelNoise Other]
            ic_classifications = EEG.etc.ic_classification.ICLabel.classifications;
            ic_labels = {'Brain', 'Muscle', 'Eye', 'Heart', 'Line Noise', 'Channel Noise', 'Other'};
            
            fprintf('IC classification completed for %d components\n', size(ic_classifications, 1));
            
        catch ME
            fprintf('ERROR: ICLabel classification failed for subject %s: %s\n', subject_id, ME.message);
            success = false;
            spec_results = [];
            return;
        end
    
        % 2. Flag and remove artifactual components
        % Initialize tracking variables for rejected ICs
        rejected_ics = [];
        rejected_ic_labels = {};
        
        if ~isempty(ic_classifications)
            fprintf('Flagging artifactual components for rejection...\n');
            
            muscle_thresh = cfg.params.ica.muscle_thresh;
            eye_thresh = cfg.params.ica.eye_thresh;
            
            try
                EEG = pop_icflag(EEG, ...
                    [NaN NaN; ...        % Brain (don't reject)
                     muscle_thresh 1; ... % Muscle (reject if > threshold)
                     eye_thresh 1; ...    % Eye (reject if > threshold)  
                     NaN NaN; ...        % Heart (don't reject)
                     NaN NaN; ...        % Line noise (don't reject)
                     NaN NaN; ...        % Channel noise (don't reject)
                     NaN NaN]);          % Other (don't reject)
                
                rejected_ics = find(EEG.reject.gcompreject);
                
                % Analyze what types of components were rejected
                if ~isempty(rejected_ics)
                    rejected_ic_labels = cell(length(rejected_ics), 1);
                    for i = 1:length(rejected_ics)
                        ic_idx = rejected_ics(i);
                        [~, max_class] = max(ic_classifications(ic_idx, :));
                        rejected_ic_labels{i} = ic_labels{max_class};
                    end
                    
                    fprintf('Flagged %d/%d components for rejection: %s\n', ...
                            length(rejected_ics), size(ic_classifications, 1), mat2str(rejected_ics));
                    
                    [unique_types, ~, idx] = unique(rejected_ic_labels);
                    type_counts = accumarray(idx, 1);
                    for i = 1:length(unique_types)
                        fprintf('  %s: %d components\n', unique_types{i}, type_counts(i));
                    end
                else
                    fprintf('No components flagged for rejection.\n');
                end
                
            catch ME
                fprintf('ERROR: Component flagging failed for subject %s: %s\n', subject_id, ME.message);
                success = false;
                spec_results = [];
                return;
            end
        else
            fprintf('No IC classification available - skipping automatic rejection.\n');
        end
        
        % 3. Remove artifactual components
        if ~isempty(rejected_ics)
            fprintf('Removing %d artifactual independent components...\n', length(rejected_ics));
            EEG = pop_subcomp(EEG, rejected_ics, 0);
            fprintf('Component removal completed.\n');
        else
            fprintf('No components flagged for removal.\n');
        end
        
        % 4. Spherical interpolation to restore 64-channel montage
        fprintf('Performing spherical interpolation to restore full montage...\n');
        
        % Capture channel info before interpolation
        n_chans_before_interp = EEG.nbchan;
        channels_before_interp = {EEG.chanlocs.labels};
        
        % Load original channel locations for interpolation reference
        if exist(cfg.channels.chan_locs, 'file')
            chanlocs = load(cfg.channels.chan_locs);
            EEG = pop_interp(EEG, chanlocs.chan_locs, 'spherical');
            fprintf('Interpolation completed: %d channels -> %d channels\n', n_chans_before_interp, EEG.nbchan);
        else
            fprintf('ERROR: Channel locations file not found for subject %s: %s\n', subject_id, cfg.channels.chan_locs);
            success = false;
            spec_results = [];
            return;
        end
        
        % Capture channel labels after interpolation
        channels_after_interp = {EEG.chanlocs.labels};
        
        % 5. Re-reference to average signal
        fprintf('Re-referencing to average signal...\n');
        EEG = pop_reref(EEG, []);
        fprintf('Average referencing completed.\n');
        
        % Save final processed EEG data
        fprintf('Saving final processed EEG data...\n');
        final_eeg_filename = sprintf('%s-final.set', subject_id);
        pop_saveset(EEG, 'filename', final_eeg_filename, 'filepath', fullfile(cfg.paths.output, 'spectral'));
        fprintf('Final EEG data saved: %s\n', fullfile(cfg.paths.output, 'spectral', final_eeg_filename));
    
        % Update matrix dimensions after interpolation
        n_chans = EEG.nbchan;
        
        % Initialize result matrices with final channel count
        spectra_db = NaN(n_chans, n_freqs, n_blocks);      % PSD in dB
        spectra_psd = NaN(n_chans, n_freqs, n_blocks);     % PSD in μV²/Hz
        freqs = NaN(n_freqs, n_blocks);                     % Frequency bins
        paf_channels = NaN(n_chans, n_blocks);              % Peak alpha frequency per channel
        cog_channels = NaN(n_chans, n_blocks);              % Center of gravity per channel
        
        % Update dataset name
        EEG.setname = sprintf('%s-final', subject_id);
        EEG = eeg_checkset(EEG);
        
        fprintf('Processing %d blocks with %d channels...\n', n_blocks, n_chans);
    
        % 6. Process each block individually
        for b = 1:n_blocks
        block_name = block_names{b};
        fprintf('\nProcessing block %d/%d: %s (%s)\n', b, n_blocks, block_name, block_descriptions{b});
        
        try
            % Epoch around this specific block
            this_epoch = [0 60];  % 60-second epochs as in original
            fprintf('  Epoching around event %s with window [%d %d] seconds...\n', block_name, this_epoch(1), this_epoch(2));
            
            EEG_block = pop_epoch(EEG, {block_name}, this_epoch, 'epochinfo', 'yes');
            
            if isempty(EEG_block.data)
                fprintf('  WARNING: Block %s resulted in empty data, skipping.\n', block_name);
                continue;
            end
            
            fprintf('  Epoched data: %d samples (%.1f seconds)\n', EEG_block.pnts, EEG_block.pnts/EEG_block.srate);
            
            % Spectral decomposition using spectopo
            fprintf('  Computing PSD using spectopo...\n');
            [spectra_db(:,:,b), freqs(:,b)] = spectopo(...
                EEG_block.data, ...
                0, ...                          % frames per epoch (0 = use all data)
                EEG_block.srate, ...           % sampling rate
                'winsize', wsize*EEG_block.srate, ... % window size in samples
                'overlap', overlap*EEG_block.srate, ... % overlap in samples
                'plot', 'off');                % don't plot
            
            % Convert from dB to PSD (μV²/Hz)
            spectra_psd(:,:,b) = 10.^(spectra_db(:,:,b)/10);
            
            fprintf('  PSD computed: %d frequency bins from %.2f to %.2f Hz\n', ...
                    length(freqs(:,b)), freqs(1,b), freqs(end,b));
            
            % Individual Alpha Frequency (IAF) analysis using restingIAF
            if exist('restingIAF', 'file')
                fprintf('  Computing Individual Alpha Frequency (IAF)...\n');
                
                try
                    % Calculate number of samples for consistent window with spectopo
                    tlen = N;  % Same window size as spectral analysis
                    
                    [pSum, pChans, f_iaf] = restingIAF(...
                        EEG_block.data, ...        % EEG data
                        EEG_block.nbchan, ...      % number of channels
                        cmin, ...                  % min channels for grand average
                        frange, ...                % spectral range
                        EEG_block.srate, ...       % sampling rate
                        alpha_window, ...          % alpha peak search window
                        fw, ...                    % SGF frame width
                        poly, ...                  % SGF polynomial order
                        'taper', 'hamming', ...    % taper type
                        'tlen', tlen, ...          % window size (samples)
                        'mdiff', 0.2);             % peak must be 20% larger
                    
                    % Store IAF results
                    if ~isempty(pChans)
                        paf_channels(:,b) = [pChans.peaks];  % Peak alpha freq per channel
                        cog_channels(:,b) = [pChans.gravs];  % Center of gravity per channel
                    end
                    
                    if ~isempty(pSum)
                        iaf_summary(b,:) = [pSum.paf pSum.cog];  % Summary IAF
                        fprintf('  IAF computed: PAF=%.2f Hz, COG=%.2f Hz\n', pSum.paf, pSum.cog);
                    else
                        fprintf('  WARNING: Could not compute summary IAF for block %s\n', block_name);
                    end
                    
                catch ME
                    fprintf('  WARNING: IAF computation failed for block %s: %s\n', block_name, ME.message);
                end
            else
                fprintf('ERROR: restingIAF function not found for subject %s\n', subject_id);
                success = false;
                spec_results = [];
                return;
            end
            
            fprintf('  Block %s processing completed.\n', block_name);
            
        catch ME
            fprintf('  ERROR: Failed to process block %s: %s\n', block_name, ME.message);
            continue;
        end
        end
    
        % 7. Store results in output structure
        spec_results.spectra_db = spectra_db;         % PSD in dB
        spec_results.spectra_psd = spectra_psd;       % PSD in μV²/Hz  
        spec_results.freqs = freqs;                   % Frequency bins
        spec_results.paf = paf_channels;              % Peak alpha frequency per channel
        spec_results.cog = cog_channels;              % Center of gravity per channel
        spec_results.iaf = iaf_summary;               % Summary IAF per block [PAF, COG]
        
        % 8. Calculate summary statistics
        n_successful_blocks = sum(~isnan(squeeze(spectra_db(1,1,:))));  % Count non-NaN blocks
        mean_psd_power = nanmean(spectra_psd(:), 'all');
        
        % Store metadata as flat fields
        spec_results.meta_original_channels = orig_n_chans;
        spec_results.meta_final_channels = n_chans;
        spec_results.meta_channels_before_interp = channels_before_interp;
        spec_results.meta_channels_after_interp = channels_after_interp;
        spec_results.meta_rejected_ics = rejected_ics;
        spec_results.meta_rejected_ic_labels = rejected_ic_labels;
        spec_results.meta_n_blocks_requested = n_blocks;
        spec_results.meta_n_blocks_processed = n_successful_blocks;
        spec_results.meta_n_frequency_bins = size(freqs, 1);
        spec_results.meta_frequency_resolution = freqs(2,1) - freqs(1,1);  % Hz per bin
        spec_results.meta_mean_psd_power = mean_psd_power;
        spec_results.meta_processing_complete = posixtime(datetime('now'));
        
        % 9. Save spectral results to output/spectral directory
        fprintf('Saving spectral analysis results for subject %s...\n', subject_id);
        spectral_filename = sprintf('%s-spectral.mat', subject_id);
        spectral_filepath = fullfile(cfg.paths.output, 'spectral', spectral_filename);
        save(spectral_filepath, 'spec_results');
        fprintf('Spectral results saved: %s\n', spectral_filepath);
        
        % Set success flag
        success = true;
        
        % Final summary
        fprintf('\nSpectral analysis completed for subject %s:\n', subject_id);
        fprintf('  Original channels: %d, Final channels: %d\n', orig_n_chans, n_chans);
        fprintf('  Blocks processed: %d/%d\n', n_successful_blocks, n_blocks);
        fprintf('  Frequency bins: %d (%.3f Hz resolution)\n', size(freqs, 1), spec_results.meta_frequency_resolution);
        fprintf('  Mean PSD power: %.2e μV²/Hz\n', mean_psd_power);
        fprintf('Spectral analysis processing completed successfully.\n');
        
    catch ME
        fprintf('ERROR: Spectral analysis failed for subject %s: %s\n', subject_id, ME.message);
        if ~isempty(ME.stack)
            fprintf('Error occurred at: %s\n', ME.stack(1).name);
        end
        success = false;
        spec_results = [];
        return;
    end

end