function EEG_blocks = extract_blocks(EEG, cfg)
% EXTRACT_BLOCKS - Extract resting state blocks from preprocessed EEG data
%
% Usage:
%   EEG_blocks = extract_blocks(EEG, cfg)
%
% Inputs:
%   EEG - EEGLAB EEG structure (preprocessed)
%   cfg - Configuration structure from config()
%
% Outputs:
%   EEG_blocks - EEGLAB EEG structure with only resting state blocks
%
% This function extracts the resting state blocks defined in cfg.blocks.names
% using the epoch window specified in cfg.blocks.epoch_window. This is done
% BEFORE ICA to avoid contamination from stimulation artifacts.

    subject_id = extractBefore(EEG.setname, '-prepro');
    fprintf('Extracting resting state blocks for subject %s...\n', subject_id);
    
    % Get block definitions
    block_names = cfg.blocks.names;
    epoch_window = cfg.blocks.epoch_window;  % e.g., [-1 61] seconds
    
    fprintf('Searching for blocks: %s\n', strjoin(block_names, ', '));
    fprintf('Using epoch window: [%d %d] seconds\n', epoch_window(1), epoch_window(2));
    
    % Check which events exist in the data
    if isempty(EEG.event)
        error('No events found in EEG data for subject %s', subject_id);
    end
    
    % Get all event types
    event_types = {EEG.event.type};
    if isnumeric([EEG.event.type])
        event_types = arrayfun(@num2str, [EEG.event.type], 'UniformOutput', false);
    end
    
    % Find which blocks are actually present in the data
    blocks_found = {};
    blocks_missing = {};
    
    for i = 1:length(block_names)
        if any(strcmp(event_types, block_names{i}))
            blocks_found{end+1} = block_names{i};
        else
            blocks_missing{end+1} = block_names{i};
        end
    end
    
    if isempty(blocks_found)
        error('None of the specified blocks found in data for subject %s. Available events: %s', ...
              subject_id, strjoin(unique(event_types), ', '));
    end
    
    fprintf('Found %d/%d blocks: %s\n', length(blocks_found), length(block_names), strjoin(blocks_found, ', '));
    if ~isempty(blocks_missing)
        fprintf('Missing blocks: %s\n', strjoin(blocks_missing, ', '));
    end
    
    % Extract the resting state blocks using pop_rmdat
    try
        fprintf('Extracting blocks using pop_rmdat...\n');
        EEG_blocks = pop_rmdat(EEG, blocks_found, epoch_window, 0);
        
        % Check if extraction was successful
        if isempty(EEG_blocks.data)
            error('Block extraction resulted in empty data');
        end
        
        fprintf('Block extraction successful:\n');
        fprintf('  Original data: %d samples (%.1f seconds)\n', EEG.pnts, EEG.pnts/EEG.srate);
        fprintf('  Extracted data: %d samples (%.1f seconds)\n', EEG_blocks.pnts, EEG_blocks.pnts/EEG_blocks.srate);
        fprintf('  Reduction: %.1f%% of original data retained\n', 100 * EEG_blocks.pnts / EEG.pnts);
        
    catch ME
        error('Block extraction failed for subject %s: %s', subject_id, ME.message);
    end
    
    % Update dataset name
    EEG_blocks.setname = sprintf('%s-blocks', subject_id);
    EEG_blocks = eeg_checkset(EEG_blocks);
    
    % Validate extracted data quality
    if any(isnan(EEG_blocks.data(:))) || any(isinf(EEG_blocks.data(:)))
        warning('Extracted data contains NaN or Inf values for subject %s', subject_id);
    end
    
    data_range = [min(EEG_blocks.data(:)), max(EEG_blocks.data(:))];
    fprintf('Extracted data range: [%.2f, %.2f] μV\n', data_range(1), data_range(2));
    
    fprintf('Block extraction completed for subject %s.\n', subject_id);
    
end