function EEG_interp = interpolate_channels(EEG_clean, bad_channels, orig_chanlocs, cfg)
% INTERPOLATE_CHANNELS - Interpolate bad channels back to restore full montage
%
% Usage:
%   EEG_interp = interpolate_channels(EEG_clean, bad_channels, orig_chanlocs, cfg)
%
% Inputs:
%   EEG_clean     - EEGLAB EEG structure after ICA cleaning (missing bad channels)
%   bad_channels  - Vector of bad channel indices from original montage
%   orig_chanlocs - Original channel locations structure (before bad channel removal)
%   cfg           - Configuration structure from config()
%
% Outputs:
%   EEG_interp - EEGLAB EEG structure with bad channels interpolated back
%
% This function restores the full electrode montage by interpolating the 
% bad channels that were removed before ICA using spherical interpolation.

    subject_id = extractBefore(EEG_clean.setname, '-ica-clean');
    fprintf('Interpolating channels for subject %s...\n', subject_id);
    
    % If no bad channels, just update the name and return
    if isempty(bad_channels)
        fprintf('No bad channels to interpolate.\n');
        EEG_interp = EEG_clean;
        EEG_interp.setname = sprintf('%s-final', subject_id);
        return;
    end
    
    fprintf('Interpolating %d bad channels: %s\n', length(bad_channels), mat2str(bad_channels));
    
    % Store current state for comparison
    pre_interp_chans = EEG_clean.nbchan;
    
    % Create a temporary EEG structure with original channel locations
    % This is needed for pop_interp to work correctly
    EEG_temp = EEG_clean;
    
    % We need to add the missing channels back as empty channels first
    orig_n_chans = length(orig_chanlocs);
    
    % Create full data matrix with NaN for missing channels
    full_data = nan(orig_n_chans, EEG_clean.pnts);
    
    % Fill in the existing (good) channels
    good_channels = setdiff(1:orig_n_chans, bad_channels);
    if length(good_channels) ~= EEG_clean.nbchan
        error('Mismatch between expected good channels (%d) and actual channels (%d)', ...
              length(good_channels), EEG_clean.nbchan);
    end
    
    full_data(good_channels, :) = EEG_clean.data;
    
    % Create temporary EEG structure with full montage
    EEG_temp.data = full_data;
    EEG_temp.nbchan = orig_n_chans;
    EEG_temp.chanlocs = orig_chanlocs;
    
    % Update channel info to match original
    if isfield(EEG_clean, 'chaninfo') && ~isempty(EEG_clean.chaninfo)
        EEG_temp.chaninfo.ndchanlocs = orig_n_chans;
    end
    
    fprintf('Created temporary structure: %d channels (%d good + %d bad)\n', ...
            EEG_temp.nbchan, length(good_channels), length(bad_channels));
    
    % Perform spherical interpolation
    try
        fprintf('Performing spherical interpolation...\n');
        EEG_interp = pop_interp(EEG_temp, bad_channels, 'spherical');
        
        fprintf('Interpolation completed successfully.\n');
        
        % Verify interpolation worked
        if EEG_interp.nbchan ~= orig_n_chans
            error('Interpolation failed: expected %d channels, got %d', orig_n_chans, EEG_interp.nbchan);
        end
        
        % Check that interpolated channels no longer contain NaN
        interp_data_check = EEG_interp.data(bad_channels, :);
        if any(isnan(interp_data_check(:)))
            warning('Some interpolated channels still contain NaN values');
        else
            fprintf('All interpolated channels contain valid data.\n');
        end
        
    catch ME
        error('Channel interpolation failed for subject %s: %s', subject_id, ME.message);
    end
    
    % Update processing history
    if ~isfield(EEG_interp.etc, 'processing_history')
        EEG_interp.etc.processing_history = EEG_clean.etc.processing_history;
    end
    hist_idx = length(EEG_interp.etc.processing_history) + 1;
    EEG_interp.etc.processing_history{hist_idx} = sprintf('Interpolated bad channels: %s', datetime('now'));
    EEG_interp.etc.processing_history{hist_idx+1} = sprintf('Channels interpolated: %s', mat2str(bad_channels));
    EEG_interp.etc.processing_history{hist_idx+2} = sprintf('Final montage: %d channels restored', EEG_interp.nbchan);
    
    % Store interpolation information
    EEG_interp.etc.interpolation_info = struct();
    EEG_interp.etc.interpolation_info.bad_channels = bad_channels;
    EEG_interp.etc.interpolation_info.pre_interp_channels = pre_interp_chans;
    EEG_interp.etc.interpolation_info.post_interp_channels = EEG_interp.nbchan;
    EEG_interp.etc.interpolation_info.method = 'spherical';
    EEG_interp.etc.interpolation_info.timestamp = datetime('now');
    
    if ~isempty(orig_chanlocs) && isfield(orig_chanlocs, 'labels')
        bad_labels = {orig_chanlocs(bad_channels).labels};
        EEG_interp.etc.interpolation_info.bad_channel_labels = bad_labels;
        fprintf('Interpolated channel labels: %s\n', strjoin(bad_labels, ', '));
    end
    
    % Update dataset name
    EEG_interp.setname = sprintf('%s-final', subject_id);
    EEG_interp = eeg_checkset(EEG_interp);
    
    % Final validation
    data_range = [min(EEG_interp.data(:)), max(EEG_interp.data(:))];
    
    % Calculate statistics for interpolated vs original channels
    interp_data_std = std(EEG_interp.data(bad_channels, :), 0, 2);
    good_data_std = std(EEG_interp.data(good_channels, :), 0, 2);
    
    fprintf('\nChannel interpolation completed for subject %s:\n', subject_id);
    fprintf('  Channels before interpolation: %d\n', pre_interp_chans);
    fprintf('  Channels after interpolation: %d\n', EEG_interp.nbchan);
    fprintf('  Bad channels interpolated: %d\n', length(bad_channels));
    fprintf('  Final data range: [%.2f, %.2f] μV\n', data_range(1), data_range(2));
    fprintf('  Interpolated channels std: %.2f ± %.2f μV\n', mean(interp_data_std), std(interp_data_std));
    fprintf('  Original channels std: %.2f ± %.2f μV\n', mean(good_data_std), std(good_data_std));
    fprintf('  Full electrode montage restored.\n');
    
end