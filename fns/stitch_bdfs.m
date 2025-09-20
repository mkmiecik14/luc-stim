function stitch_bdfs(subject_id, part1_filename, part2_filename, outname, cfg)
% STITCH_BDFS - Combine two BDF files for a subject into one continuous recording
%
% Usage:
%   stitch_bdfs('1461831842001', '1461831842001_part1.bdf', '1461831842001_part2.bdf', '1461831842001_1', cfg)
%
% This function loads two BDF files, concatenates them, and saves as a single BDF
% that can be processed by the normal pipeline
    
    % Validate that outname exists and is a string/char
    validateattributes(outname, {'char', 'string'}, {'nonempty'});

    fprintf('Stitching BDF files for subject %s...\n', subject_id);
    
    % Load first BDF file
    bdf1_path = fullfile(cfg.paths.data, part1_filename);
    if ~exist(bdf1_path, 'file')
        error('First BDF file not found: %s', bdf1_path);
    end
    
    fprintf('Loading first BDF: %s\n', part1_filename);
    EEG1 = pop_biosig(bdf1_path, ...
        'ref', [1], ...
        'refoptions', {'keepref','on'}, ...
        'importannot', 'off', ...
        'bdfeventmode', 6);
    
    % Load second BDF file
    bdf2_path = fullfile(cfg.paths.data, part2_filename);
    if ~exist(bdf2_path, 'file')
        error('Second BDF file not found: %s', bdf2_path);
    end
    
    fprintf('Loading second BDF: %s\n', part2_filename);
    EEG2 = pop_biosig(bdf2_path, ...
        'ref', [1], ...
        'refoptions', {'keepref','on'}, ...
        'importannot', 'off', ...
        'bdfeventmode', 6);
    
    % Check that files have same number of channels and sampling rate
    if EEG1.nbchan ~= EEG2.nbchan
        error('Channel count mismatch: File 1 has %d channels, File 2 has %d channels', ...
              EEG1.nbchan, EEG2.nbchan);
    end
    
    if EEG1.srate ~= EEG2.srate
        error('Sampling rate mismatch: File 1 is %.1f Hz, File 2 is %.1f Hz', ...
              EEG1.srate, EEG2.srate);
    end
    
    fprintf('Both files validated: %d channels, %.1f Hz\n', EEG1.nbchan, EEG1.srate);
    
    % Merge the datasets using EEGLAB function
    fprintf('Merging datasets...\n');
    EEG_combined = pop_mergeset(EEG1, EEG2, 1);
    
    % Update dataset info
    EEG_combined.setname = sprintf('%s-stitched', subject_id);
    EEG_combined = eeg_checkset(EEG_combined);
    
    % Save the combined dataset as a .set file (EEGLAB format)
    % This allows your pipeline to load it normally
    output_filename = sprintf('%s.set', outname);
    output_path = fullfile(cfg.paths.output, 'raw-data');
    
    fprintf('Saving stitched dataset: %s\n', output_filename);
    EEG_combined = pop_saveset(EEG_combined, 'filename', output_filename, 'filepath', output_path);
    
    fprintf('Successfully stitched and saved: %s\n', output_path);
    fprintf('Original file 1: %d samples (%.1f seconds)\n', EEG1.pnts, EEG1.pnts/EEG1.srate);
    fprintf('Original file 2: %d samples (%.1f seconds)\n', EEG2.pnts, EEG2.pnts/EEG2.srate);
    fprintf('Combined file: %d samples (%.1f seconds)\n', EEG_combined.pnts, EEG_combined.pnts/EEG_combined.srate);
    
end