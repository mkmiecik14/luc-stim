function spec_results = spectral_analysis(EEG_final, cfg)
% SPECTRAL_ANALYSIS - Perform spectral decomposition and IAF analysis
%
% Usage:
%   spec_results = spectral_analysis(EEG_final, cfg)
%
% Inputs:
%   EEG_final - EEGLAB EEG structure (final processed data with full montage)
%   cfg       - Configuration structure from config()
%
% Outputs:
%   spec_results - Structure containing spectral analysis results
%
% This function performs spectral analysis on individual resting state blocks:
% 1. Epoch data by blocks (eyes open/closed conditions)
% 2. Calculate power spectral density (PSD) using spectopo
% 3. Compute Individual Alpha Frequency (IAF) using restingIAF
% 4. Store results in organized structure
%
% Note: Pink/white noise correction is NOT applied in this simplified version

    subject_id = extractBefore(EEG_final.setname, '-final');
    fprintf('Performing spectral analysis for subject %s...\n', subject_id);
    
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
    spec_results.timestamp = datetime('now');
    spec_results.parameters = struct('wsize', wsize, 'overlap', overlap, 'frange', frange, ...
                                   'alpha_window', alpha_window, 'cmin', cmin, 'fw', fw, 'poly', poly);
    spec_results.block_names = block_names;
    spec_results.block_descriptions = block_descriptions;
    
    % Pre-allocate result matrices
    n_blocks = length(block_names);
    n_chans = EEG_final.nbchan;
    
    % Calculate expected frequency resolution
    N = EEG_final.srate * wsize;                    % Number of data points in FFT window
    n_freqs = floor(N/2) + 1;                       % Number of frequency bins
    
    % Initialize matrices
    spectra_db = NaN(n_chans, n_freqs, n_blocks);      % PSD in dB
    spectra_psd = NaN(n_chans, n_freqs, n_blocks);     % PSD in μV²/Hz
    freqs = NaN(n_freqs, n_blocks);                     % Frequency bins
    paf_channels = NaN(n_chans, n_blocks);              % Peak alpha frequency per channel
    cog_channels = NaN(n_chans, n_blocks);              % Center of gravity per channel  
    iaf_summary = NaN(n_blocks, 2);                     % Summary IAF [PAF, COG] per block

    % 4. Classify independent components
        fprintf('Classifying independent components using ICLabel...\n');
        
        try
            EEG = pop_iclabel(EEG, 'default');
            
            % ICLabel classifications: [Brain Muscle Eye Heart LineNoise ChannelNoise Other]
            ic_classifications = EEG.etc.ic_classification.ICLabel.classifications;
            ic_labels = {'Brain', 'Muscle', 'Eye', 'Heart', 'Line Noise', 'Channel Noise', 'Other'};
            
            fprintf('IC classification completed for %d components\n', size(ic_classifications, 1));
            
        catch ME
            warning('ICLabel classification failed: %s. Continuing without IC classification.', ME.message);
            ic_classifications = [];
        end
    
        % 5. Flag and remove artifactual components
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
                    rejected_types = cell(length(rejected_ics), 1);
                    for i = 1:length(rejected_ics)
                        ic_idx = rejected_ics(i);
                        [~, max_class] = max(ic_classifications(ic_idx, :));
                        rejected_types{i} = ic_labels{max_class};
                    end
                    
                    fprintf('Flagged %d/%d components for rejection: %s\n', ...
                            length(rejected_ics), size(ic_classifications, 1), mat2str(rejected_ics));
                    
                    [unique_types, ~, idx] = unique(rejected_types);
                    type_counts = accumarray(idx, 1);
                    for i = 1:length(unique_types)
                        fprintf('  %s: %d components\n', unique_types{i}, type_counts(i));
                    end
                else
                    fprintf('No components flagged for rejection.\n');
                end
                
            catch ME
                warning(sprintf('Component flagging failed: %s', ME.message));
                rejected_ics = [];
            end
        else
            rejected_ics = [];
            fprintf('No IC classification available - skipping automatic rejection.\n');
        end
        
        % 6. Remove artifactual components
        if ~isempty(rejected_ics)
            fprintf('Removing %d artifactual independent components...\n', length(rejected_ics));
            EEG = pop_subcomp(EEG, rejected_ics, 0);
            fprintf('Component removal completed.\n');
        else
            fprintf('No components flagged for removal.\n');
        end
    
    fprintf('Processing %d blocks with %d channels...\n', n_blocks, n_chans);
    
    % Process each block individually
    for b = 1:n_blocks
        block_name = block_names{b};
        fprintf('\nProcessing block %d/%d: %s (%s)\n', b, n_blocks, block_name, block_descriptions{b});
        
        try
            % Epoch around this specific block
            this_epoch = [0 60];  % 60-second epochs as in original
            fprintf('  Epoching around event %s with window [%d %d] seconds...\n', block_name, this_epoch(1), this_epoch(2));
            
            EEG_block = pop_epoch(EEG_final, {block_name}, this_epoch, 'epochinfo', 'yes');
            
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
                fprintf('  WARNING: restingIAF function not found, skipping IAF analysis.\n');
            end
            
            fprintf('  Block %s processing completed.\n', block_name);
            
        catch ME
            fprintf('  ERROR: Failed to process block %s: %s\n', block_name, ME.message);
            continue;
        end
    end
    
    % Store results in output structure
    spec_results.spectra_db = spectra_db;         % PSD in dB
    spec_results.spectra_psd = spectra_psd;       % PSD in μV²/Hz  
    spec_results.freqs = freqs;                   % Frequency bins
    spec_results.paf = paf_channels;              % Peak alpha frequency per channel
    spec_results.cog = cog_channels;              % Center of gravity per channel
    spec_results.iaf = iaf_summary;               % Summary IAF per block [PAF, COG]
    
    % Calculate summary statistics
    n_successful_blocks = sum(~isnan(squeeze(spectra_db(1,1,:))));  % Count non-NaN blocks
    mean_psd_power = nanmean(spectra_psd(:), 'all');
    
    % Store metadata
    spec_results.metadata = struct();
    spec_results.metadata.n_channels = n_chans;
    spec_results.metadata.n_blocks_requested = n_blocks;
    spec_results.metadata.n_blocks_processed = n_successful_blocks;
    spec_results.metadata.n_frequency_bins = size(freqs, 1);
    spec_results.metadata.frequency_resolution = freqs(2,1) - freqs(1,1);  % Hz per bin
    spec_results.metadata.mean_psd_power = mean_psd_power;
    spec_results.metadata.processing_complete = datetime('now');
    
    % Final summary
    fprintf('\nSpectral analysis completed for subject %s:\n', subject_id);
    fprintf('  Blocks processed: %d/%d\n', n_successful_blocks, n_blocks);
    fprintf('  Frequency bins: %d (%.3f Hz resolution)\n', size(freqs, 1), spec_results.metadata.frequency_resolution);
    fprintf('  Mean PSD power: %.2e μV²/Hz\n', mean_psd_power);
    
    if exist('restingIAF', 'file')
        valid_iaf = sum(~isnan(iaf_summary(:,1)));
        if valid_iaf > 0
            mean_paf = nanmean(iaf_summary(:,1));
            mean_cog = nanmean(iaf_summary(:,2));
            fprintf('  IAF summary: PAF=%.2f±%.2f Hz, COG=%.2f±%.2f Hz (%d blocks)\n', ...
                    mean_paf, nanstd(iaf_summary(:,1)), mean_cog, nanstd(iaf_summary(:,2)), valid_iaf);
        else
            fprintf('  No valid IAF estimates obtained.\n');
        end
    end
    
    fprintf('Spectral analysis completed.\n');
    
end