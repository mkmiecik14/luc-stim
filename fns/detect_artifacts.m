function [bad_channels, artifact_info] = detect_artifacts(EEG, cfg)
% DETECT_ARTIFACTS - Automatic detection of bad channels and artifacts
%
% Usage:
%   [bad_channels, artifact_info] = detect_artifacts(EEG, cfg)
%
% Inputs:
%   EEG - EEGLAB EEG structure (with resting state blocks extracted)
%   cfg - Configuration structure from config()
%
% Outputs:
%   bad_channels  - Vector of bad channel indices
%   artifact_info - Structure containing detailed artifact detection info
%
% This function uses multiple criteria to automatically detect bad channels:
% 1. High variance channels (statistical outliers)
% 2. Low correlation with neighboring channels
% 3. High-frequency noise indicators
% 4. Flat or near-flat channels

    subject_id = extractBefore(EEG.setname, '-blocks');
    fprintf('Detecting artifacts for subject %s...\n', subject_id);
    
    % Initialize outputs
    bad_channels = [];
    artifact_info = struct();
    
    % Get parameters
    var_thresh = cfg.params.bad_chan.var_thresh;      % Standard deviations for variance threshold
    corr_thresh = cfg.params.bad_chan.corr_thresh;    % Minimum correlation with neighbors
    hf_thresh = cfg.params.bad_chan.hf_thresh;        % High-frequency noise threshold
    
    n_chans = EEG.nbchan;
    n_samples = EEG.pnts;
    
    fprintf('Analyzing %d channels with %d samples...\n', n_chans, n_samples);
    
    % 1. VARIANCE-BASED DETECTION
    fprintf('1. Detecting high-variance channels...\n');
    
    % Calculate variance for each channel
    chan_var = var(EEG.data, 0, 2);  % Variance across time for each channel
    
    % Find outliers using robust statistics
    median_var = median(chan_var);
    mad_var = mad(chan_var, 1);  % Median absolute deviation
    var_threshold = median_var + var_thresh * mad_var;
    
    bad_var = find(chan_var > var_threshold);
    
    fprintf('   Median variance: %.2f, MAD: %.2f, Threshold: %.2f\n', median_var, mad_var, var_threshold);
    fprintf('   High-variance channels (%d): %s\n', length(bad_var), mat2str(bad_var));
    
    % 2. CORRELATION-BASED DETECTION
    fprintf('2. Detecting low-correlation channels...\n');
    
    % Calculate correlation matrix
    corr_matrix = corrcoef(EEG.data');
    
    % For each channel, calculate mean correlation with all other channels
    mean_corr = zeros(n_chans, 1);
    for i = 1:n_chans
        % Exclude self-correlation (diagonal)
        other_chans = setdiff(1:n_chans, i);
        mean_corr(i) = mean(abs(corr_matrix(i, other_chans)), 'omitnan');
    end
    
    bad_corr = find(mean_corr < corr_thresh);
    
    fprintf('   Mean correlation threshold: %.2f\n', corr_thresh);
    fprintf('   Low-correlation channels (%d): %s\n', length(bad_corr), mat2str(bad_corr));
    
    % 3. HIGH-FREQUENCY NOISE DETECTION
    fprintf('3. Detecting high-frequency noise...\n');
    
    % Calculate power in high-frequency range (> 50 Hz)
    hf_start_freq = 50;  % Hz
    if EEG.srate/2 > hf_start_freq
        % Use FFT to get power in high frequencies
        n_fft = min(2^nextpow2(n_samples), 2^14);  % Limit FFT size for efficiency
        freqs = (0:n_fft/2) * EEG.srate / n_fft;
        hf_idx = freqs > hf_start_freq;
        
        hf_power = zeros(n_chans, 1);
        for i = 1:n_chans
            fft_data = fft(EEG.data(i, :), n_fft);
            power_spectrum = abs(fft_data(1:n_fft/2+1)).^2;
            hf_power(i) = mean(power_spectrum(hf_idx));
        end
        
        % Find channels with excessive high-frequency power
        median_hf = median(hf_power);
        mad_hf = mad(hf_power, 1);
        hf_threshold = median_hf + hf_thresh * mad_hf;
        
        bad_hf = find(hf_power > hf_threshold);
        
        fprintf('   Median HF power: %.2e, MAD: %.2e, Threshold: %.2e\n', median_hf, mad_hf, hf_threshold);
        fprintf('   High-frequency noise channels (%d): %s\n', length(bad_hf), mat2str(bad_hf));
    else
        bad_hf = [];
        fprintf('   Sampling rate too low for HF noise detection (Nyquist: %.1f Hz)\n', EEG.srate/2);
    end
    
    % 4. FLAT CHANNEL DETECTION
    fprintf('4. Detecting flat channels...\n');
    
    % Calculate range (max - min) for each channel
    chan_range = max(EEG.data, [], 2) - min(EEG.data, [], 2);
    
    % Channels with very small range are likely flat
    flat_threshold = 1.0;  % μV, very conservative threshold
    bad_flat = find(chan_range < flat_threshold);
    
    fprintf('   Range threshold: %.1f μV\n', flat_threshold);
    fprintf('   Flat channels (%d): %s\n', length(bad_flat), mat2str(bad_flat));
    
    % 5. COMBINE ALL CRITERIA
    fprintf('5. Combining all detection criteria...\n');
    
    % Collect all bad channels (union of all methods)
    all_bad = unique([bad_var; bad_corr; bad_hf; bad_flat]);
    
    % Store detailed information
    artifact_info.subject_id = subject_id;
    artifact_info.n_channels = n_chans;
    artifact_info.n_samples = n_samples;
    artifact_info.bad_variance = bad_var;
    artifact_info.bad_correlation = bad_corr;
    artifact_info.bad_hf_noise = bad_hf;
    artifact_info.bad_flat = bad_flat;
    artifact_info.variance_stats = struct('median', median_var, 'mad', mad_var, 'threshold', var_threshold);
    artifact_info.correlation_stats = struct('threshold', corr_thresh, 'mean_correlations', mean_corr);
    if exist('hf_power', 'var')
        artifact_info.hf_stats = struct('median', median_hf, 'mad', mad_hf, 'threshold', hf_threshold, 'power', hf_power);
    end
    artifact_info.flat_stats = struct('threshold', flat_threshold, 'ranges', chan_range);
    artifact_info.detection_params = cfg.params.bad_chan;
    artifact_info.timestamp = datetime('now');
    
    % Return bad channels
    bad_channels = all_bad;
    
    % Summary
    fprintf('\nArtifact detection summary for subject %s:\n', subject_id);
    fprintf('  High variance: %d channels\n', length(bad_var));
    fprintf('  Low correlation: %d channels\n', length(bad_corr));
    fprintf('  High-frequency noise: %d channels\n', length(bad_hf));
    fprintf('  Flat channels: %d channels\n', length(bad_flat));
    fprintf('  TOTAL BAD CHANNELS: %d/%d (%.1f%%)\n', length(bad_channels), n_chans, 100*length(bad_channels)/n_chans);
    
    if ~isempty(bad_channels)
        fprintf('  Bad channel indices: %s\n', mat2str(bad_channels));
        if ~isempty(EEG.chanlocs)
            bad_labels = {EEG.chanlocs(bad_channels).labels};
            fprintf('  Bad channel labels: %s\n', strjoin(bad_labels, ', '));
            artifact_info.bad_channel_labels = bad_labels;
        end
    else
        fprintf('  No bad channels detected.\n');
    end
    
    fprintf('Artifact detection completed for subject %s.\n', subject_id);
    
end