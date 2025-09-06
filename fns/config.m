function cfg = config()
% CONFIG - Initialize EEG processing environment and parameters
%
% This function replaces workspace_prep.m with a more robust, 
% machine-independent configuration system.
%
% Usage:
%   cfg = config();
%
% Environment Variables:
%   LUC_STIM_DATA_PATH - Optional custom data directory path
%                        If not set, defaults to './data' relative to project
%
% Outputs:
%   cfg - Configuration structure containing all paths, parameters, and data
%
% Matt Kmiecik
% Updated for modular pipeline

    % Get current working directory (machine-independent)
    current_dir = pwd;
    
    % Navigate to project root (assuming we're in src/ subdirectory)
    if contains(current_dir, 'src')
        main_dir = fileparts(current_dir);
    else
        main_dir = current_dir;
    end
    
    % Set working directory
    cd(main_dir);
    
    % Directory paths
    cfg.paths.main = main_dir;
    cfg.paths.fns = fullfile(main_dir, 'fns');
    cfg.paths.doc = fullfile(main_dir, 'doc');
    cfg.paths.output = fullfile(main_dir, 'output');
    
    % Data directory - check for custom path via environment variable
    data_path_env = getenv('LUC_STIM_DATA_PATH');
    if ~isempty(data_path_env)
        if exist(data_path_env, 'dir')
            cfg.paths.data = data_path_env;
            fprintf('Using custom data directory: %s\n', data_path_env);
        else
            warning('Custom data path does not exist: %s\nFalling back to default', data_path_env);
            cfg.paths.data = fullfile(main_dir, 'data');
        end
    else
        cfg.paths.data = fullfile(main_dir, 'data');
    end
    
    % Create output subdirectories if they don't exist
    output_subdirs = {'preprocessed', 'ica', 'spectral', 'logs', 'qc', 'amica'};
    for i = 1:length(output_subdirs)
        subdir_path = fullfile(cfg.paths.output, output_subdirs{i});
        if ~exist(subdir_path, 'dir')
            mkdir(subdir_path);
            fprintf('Created directory: %s\n', subdir_path);
        end
    end
    
    % Processing parameters
    cfg.params.srate_target = 256;           % Target sampling rate (Hz)
    cfg.params.linked_mast_ref = [24 61];    % Linked mastoid reference
    cfg.params.highpass_freq = 1;            % High-pass filter frequency (Hz)
    cfg.params.line_freq = [60 120 180];     % Line noise frequency (Hz)
    cfg.params.wsize = 4;                    % FFT window size (seconds)
    cfg.params.overlap = 2;                  % FFT overlap (seconds)
    
    % Bad channel detection parameters (EEGLAB pop_clean_rawdata)
    cfg.params.bad_chan.FlatlineCriterion = 5;    % Channel: flat channels
    cfg.params.bad_chan.ChannelCriterion = 0.8;   % Channel: correlation-based
    cfg.params.bad_chan.LineNoiseCriterion = 4;   % Channel: line noise
    
    % ICA parameters
    cfg.params.ica.muscle_thresh = 0.8;      % Muscle IC rejection threshold
    cfg.params.ica.eye_thresh = 0.8;         % Eye IC rejection threshold
    cfg.params.ica.extended = 1;             % Extended ICA
    cfg.params.ica.amica_path = fullfile(main_dir, 'output', 'amica');
    
    % Spectral analysis parameters
    cfg.params.spectral.frange = [1 40];     % Frequency range for analysis
    cfg.params.spectral.alpha_window = [6 15]; % Alpha peak search window
    cfg.params.spectral.cmin = 3;            % Min channels for cross-channel avg
    cfg.params.spectral.fw = 23;             % SGF frame width
    cfg.params.spectral.poly = 5;            % SGF polynomial order
    
    % Resting state block definitions
    cfg.blocks.names = {'103','105','107', '109', '203', '205', '207', '209'};
    cfg.blocks.descriptions = { ...
        'eyes open 183_PRE_RS', ... 
        'eyes closed 183_PRE_RS', ...
        'eyes open 183_PRE_RS', ...
        'eyes closed 183_PRE_RS', ...
        'eyes open 184_POST_RS', ...
        'eyes closed 184_POST_RS',  ...
        'eyes open 184_POST_RS', ...
        'eyes closed 184_POST_RS' ...
    };
    cfg.blocks.epoch_window = [-1 61];       % Epoch window around trigger (sec)
    
    % Channel setup
    cfg.channels.chan_info = fullfile(cfg.paths.doc, 'chan_info_nose_along_fixed.mat');
    cfg.channels.chan_locs = fullfile(cfg.paths.doc, 'chan_locs_nose_along_fixed.mat');
    cfg.channels.external_remove = {'EXG1', 'EXG2', 'EXG3','EXG4','EXG5','EXG6','EXG7','EXG8'};
    
    % Generate A1-A32, B1-B32 channel names for 64-channel setup
    N = 32;
    A_chans = cell(1, N);
    B_chans = A_chans;
    for j = 1:N
        A_chans{j} = sprintf('A%d', j);
        B_chans{j} = sprintf('B%d', j);
    end
    cfg.channels.keep = [A_chans B_chans];
    
    % Initialize EEGLAB
    fprintf('Initializing EEGLAB...\n');
    [ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab('nogui');
    
    % Store EEGLAB variables in config
    cfg.eeglab.ALLEEG = ALLEEG;
    cfg.eeglab.EEG = EEG;
    cfg.eeglab.CURRENTSET = CURRENTSET;
    cfg.eeglab.ALLCOM = ALLCOM;
    
    % Load participant information
    ss_info_path = fullfile(cfg.paths.doc, 'ss-info.xlsx');
    if exist(ss_info_path, 'file')
        try
            [NUM, TXT, RAW] = xlsread(ss_info_path);
            cfg.participants.NUM = NUM;
            cfg.participants.TXT = TXT;
            cfg.participants.RAW = RAW;
            cfg.participants.ss_list = string({RAW{2:size(RAW,1),1}});
            fprintf('Loaded participant information: %d participants\n', length(cfg.participants.ss_list));
        catch ME
            warning(sprintf('Could not load participant information: %s', ME.message));
            cfg.participants = [];
        end
    else
        warning('Participant info file not found: %s', ss_info_path);
        cfg.participants = [];
    end
    
    % Timestamp
    cfg.timestamp = datetime('now');
    
    fprintf('Configuration completed successfully.\n');
    fprintf('Main directory: %s\n', cfg.paths.main);
    fprintf('Data directory: %s\n', cfg.paths.data);
    fprintf('Output directory: %s\n', cfg.paths.output);
    
end