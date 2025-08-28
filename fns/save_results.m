function save_results(EEG_final, spec_results, artifact_info, ica_info, cfg)
% SAVE_RESULTS - Save processed EEG data and analysis results
%
% Usage:
%   save_results(EEG_final, spec_results, artifact_info, ica_info, cfg)
%
% Inputs:
%   EEG_final    - Final processed EEGLAB EEG structure
%   spec_results - Spectral analysis results structure
%   artifact_info- Artifact detection information structure
%   ica_info     - ICA processing information structure  
%   cfg          - Configuration structure from config()
%
% This function saves:
% 1. Final processed EEG dataset (.set/.fdt files)
% 2. Spectral results (.mat file, compatible with original format)
% 3. Processing log with all parameters and steps
% 4. Quality control report

    subject_id = spec_results.subject_id;
    fprintf('Saving results for subject %s...\n', subject_id);
    
    % Create timestamped processing ID
    timestamp_str = datestr(now, 'yyyymmdd_HHMMSS');
    
    %% 1. SAVE FINAL EEG DATASET
    fprintf('1. Saving final EEG dataset...\n');
    
    try
        % Save to preprocessed directory with descriptive name
        eeg_filename = sprintf('%s-final-processed.set', subject_id);
        eeg_filepath = fullfile(cfg.paths.output, 'preprocessed');
        
        EEG_final = pop_saveset(EEG_final, 'filename', eeg_filename, 'filepath', eeg_filepath);
        fprintf('   EEG dataset saved: %s\n', fullfile(eeg_filepath, eeg_filename));
        
    catch ME
        warning('Failed to save EEG dataset: %s', ME.message);
    end
    
    %% 2. SAVE SPECTRAL RESULTS (Compatible with original format)
    fprintf('2. Saving spectral analysis results...\n');
    
    try
        % Create spec_res structure compatible with original pipeline
        spec_res = struct();
        spec_res.spectra = spec_results.spectra_db;      % PSD in dB (original format)
        spec_res.psd = spec_results.spectra_psd;         % PSD in μV²/Hz  
        spec_res.freqs = spec_results.freqs;             % Frequency bins
        spec_res.paf = spec_results.paf;                 % Peak alpha frequency per channel
        spec_res.cog = spec_results.cog;                 % Center of gravity per channel
        spec_res.iaf = spec_results.iaf;                 % Summary IAF per block
        
        % Additional metadata for new pipeline
        spec_res.metadata = spec_results.metadata;
        spec_res.parameters = spec_results.parameters;
        spec_res.block_names = spec_results.block_names;
        spec_res.processing_timestamp = spec_results.timestamp;
        
        % Save in spectral directory
        spec_filename = sprintf('%s-spec-res.mat', subject_id);
        spec_filepath = fullfile(cfg.paths.output, 'spectral');
        
        save(fullfile(spec_filepath, spec_filename), 'spec_res');
        fprintf('   Spectral results saved: %s\n', fullfile(spec_filepath, spec_filename));
        
    catch ME
        warning('Failed to save spectral results: %s', ME.message);
    end
    
    %% 3. SAVE COMPREHENSIVE PROCESSING LOG
    fprintf('3. Saving processing log...\n');
    
    try
        % Create comprehensive processing log
        log_data = struct();
        log_data.subject_id = subject_id;
        log_data.processing_timestamp = timestamp_str;
        log_data.pipeline_version = 'modular_v1.0';
        
        % Configuration used
        log_data.config = cfg;
        
        % Processing steps and results
        log_data.artifact_detection = artifact_info;
        log_data.ica_processing = ica_info;
        log_data.spectral_analysis = spec_results;
        
        % EEG processing history
        if isfield(EEG_final.etc, 'processing_history')
            log_data.processing_history = EEG_final.etc.processing_history;
        end
        
        % Data quality metrics
        log_data.quality_metrics = struct();
        log_data.quality_metrics.final_channels = EEG_final.nbchan;
        log_data.quality_metrics.final_samples = EEG_final.pnts;
        log_data.quality_metrics.final_duration_sec = EEG_final.pnts / EEG_final.srate;
        log_data.quality_metrics.data_range = [min(EEG_final.data(:)), max(EEG_final.data(:))];
        log_data.quality_metrics.mean_power = mean(var(EEG_final.data, 0, 2));
        
        if ~isempty(artifact_info)
            log_data.quality_metrics.bad_channels_detected = length(artifact_info.subject_id);
        end
        
        if ~isempty(ica_info)
            log_data.quality_metrics.ica_components = ica_info.n_components;
            log_data.quality_metrics.rejected_components = length(ica_info.rejected_ics);
        end
        
        % Save processing log
        log_filename = sprintf('%s-processing-log-%s.mat', subject_id, timestamp_str);
        log_filepath = fullfile(cfg.paths.output, 'logs');
        
        save(fullfile(log_filepath, log_filename), 'log_data');
        fprintf('   Processing log saved: %s\n', fullfile(log_filepath, log_filename));
        
    catch ME
        warning('Failed to save processing log: %s', ME.message);
    end
    
    %% 4. SAVE QUALITY CONTROL REPORT (Text Summary)
    fprintf('4. Saving quality control report...\n');
    
    try
        qc_filename = sprintf('%s-qc-report-%s.txt', subject_id, timestamp_str);
        qc_filepath = fullfile(cfg.paths.output, 'qc', qc_filename);
        
        % Open file for writing
        fid = fopen(qc_filepath, 'w');
        if fid == -1
            error('Cannot open QC report file for writing');
        end
        
        % Write QC report
        fprintf(fid, 'EEG PROCESSING QUALITY CONTROL REPORT\n');
        fprintf(fid, '=====================================\n\n');
        fprintf(fid, 'Subject ID: %s\n', subject_id);
        fprintf(fid, 'Processing Date: %s\n', timestamp_str);
        fprintf(fid, 'Pipeline Version: modular_v1.0\n\n');
        
        % Processing summary
        fprintf(fid, 'PROCESSING SUMMARY\n');
        fprintf(fid, '------------------\n');
        fprintf(fid, 'Final Channels: %d\n', EEG_final.nbchan);
        fprintf(fid, 'Final Samples: %d\n', EEG_final.pnts);
        fprintf(fid, 'Final Duration: %.1f seconds\n', EEG_final.pnts / EEG_final.srate);
        fprintf(fid, 'Sampling Rate: %.1f Hz\n', EEG_final.srate);
        fprintf(fid, 'Data Range: [%.2f, %.2f] μV\n', min(EEG_final.data(:)), max(EEG_final.data(:)));
        
        % Artifact detection summary
        if ~isempty(artifact_info)
            fprintf(fid, '\nARTIFACT DETECTION\n');
            fprintf(fid, '------------------\n');
            fprintf(fid, 'Bad Channels Detected: %d/%d (%.1f%%)\n', ...
                    length(artifact_info.bad_variance), artifact_info.n_channels, ...
                    100 * length(artifact_info.bad_variance) / artifact_info.n_channels);
            
            if isfield(artifact_info, 'bad_channel_labels') && ~isempty(artifact_info.bad_channel_labels)
                fprintf(fid, 'Bad Channel Labels: %s\n', strjoin(artifact_info.bad_channel_labels, ', '));
            end
            
            fprintf(fid, 'Detection Criteria:\n');
            fprintf(fid, '  High variance: %d channels\n', length(artifact_info.bad_variance));
            fprintf(fid, '  Low correlation: %d channels\n', length(artifact_info.bad_correlation));
            fprintf(fid, '  High-frequency noise: %d channels\n', length(artifact_info.bad_hf_noise));
            fprintf(fid, '  Flat channels: %d channels\n', length(artifact_info.bad_flat));
        end
        
        % ICA processing summary  
        if ~isempty(ica_info)
            fprintf(fid, '\nICA PROCESSING\n');
            fprintf(fid, '--------------\n');
            fprintf(fid, 'ICA Components: %d\n', ica_info.n_components);
            fprintf(fid, 'ICA Rank: %d\n', ica_info.ica_rank);
            fprintf(fid, 'Rejected Components: %d\n', length(ica_info.rejected_ics));
            
            if ~isempty(ica_info.rejected_ics)
                fprintf(fid, 'Rejected IC indices: %s\n', mat2str(ica_info.rejected_ics));
                if isfield(ica_info, 'rejected_types')
                    unique_types = unique(ica_info.rejected_types);
                    for i = 1:length(unique_types)
                        count = sum(strcmp(ica_info.rejected_types, unique_types{i}));
                        fprintf(fid, '  %s components: %d\n', unique_types{i}, count);
                    end
                end
            end
        end
        
        % Spectral analysis summary
        fprintf(fid, '\nSPECTRAL ANALYSIS\n');
        fprintf(fid, '-----------------\n');
        fprintf(fid, 'Blocks Processed: %d/%d\n', spec_results.metadata.n_blocks_processed, spec_results.metadata.n_blocks_requested);
        fprintf(fid, 'Frequency Bins: %d\n', spec_results.metadata.n_frequency_bins);
        fprintf(fid, 'Frequency Resolution: %.3f Hz\n', spec_results.metadata.frequency_resolution);
        fprintf(fid, 'Mean PSD Power: %.2e μV²/Hz\n', spec_results.metadata.mean_psd_power);
        
        % IAF summary if available
        if ~isempty(spec_results.iaf) && any(~isnan(spec_results.iaf(:)))
            valid_iaf = ~isnan(spec_results.iaf(:,1));
            if any(valid_iaf)
                mean_paf = mean(spec_results.iaf(valid_iaf,1));
                mean_cog = mean(spec_results.iaf(valid_iaf,2));
                fprintf(fid, 'Individual Alpha Frequency:\n');
                fprintf(fid, '  Mean PAF: %.2f Hz\n', mean_paf);
                fprintf(fid, '  Mean COG: %.2f Hz\n', mean_cog);
                fprintf(fid, '  Valid blocks: %d/%d\n', sum(valid_iaf), length(valid_iaf));
            end
        end
        
        % Processing parameters
        fprintf(fid, '\nPROCESSING PARAMETERS\n');
        fprintf(fid, '---------------------\n');
        fprintf(fid, 'High-pass filter: %.1f Hz\n', cfg.params.highpass_freq);
        fprintf(fid, 'Line noise removal: %d Hz\n', cfg.params.line_freq);
        fprintf(fid, 'Sampling rate: %d Hz\n', cfg.params.srate_target);
        fprintf(fid, 'Spectral window: %.1f s\n', cfg.params.wsize);
        fprintf(fid, 'Spectral overlap: %.1f s\n', cfg.params.overlap);
        
        fclose(fid);
        fprintf('   QC report saved: %s\n', qc_filepath);
        
    catch ME
        if exist('fid', 'var') && fid ~= -1
            fclose(fid);
        end
        warning('Failed to save QC report: %s', ME.message);
    end
    
    %% 5. SUMMARY
    fprintf('\nResults successfully saved for subject %s:\n', subject_id);
    fprintf('  EEG dataset: %s/preprocessed/\n', cfg.paths.output);
    fprintf('  Spectral results: %s/spectral/\n', cfg.paths.output);
    fprintf('  Processing log: %s/logs/\n', cfg.paths.output);
    fprintf('  QC report: %s/qc/\n', cfg.paths.output);
    fprintf('All outputs compatible with original pipeline format.\n');
    
end