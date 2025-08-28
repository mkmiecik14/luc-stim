% EEG Preprocessing Pipeline Step 1
% Matt Kmiecik
% Started 03 June 2023

workspace_prep % Prepares workspace (see src/)

% Initializes subjects for batch processing (if applicable)
ss = string({RAW{2:size(RAW,1),1}});

i=1; % for testing purposes

% Preprocessing ----
for i = 1:length(ss)

    % Creating variables ----
    this_ss = ss{i};
    this_ss_path = dir(fullfile(data_dir, strcat(this_ss, '.bdf')));
    this_ss_name = this_ss_path.name;

    % Loads in raw data using biosemi ----
    EEG = pop_biosig(...
        fullfile(this_ss_path.folder, this_ss_name),...
        'ref',[1] ,...
        'refoptions',{'keepref','on'},...
        'importannot', 'off',... % does not import EDF annotations
        'bdfeventmode', 6 ... % this event mode syncs nicely with EMSE events
        );
    
    % Checks to see if an outside event file exists 
    % (i.e., if sections of the EEG were rejected due to noise)
    %this_events = strcat(num2str(this_ss), '-events.csv');
    %if isfile(fullfile(data_dir, this_events))
     % File exists...load in the visually inspected and rejected file
     % Loads in raw data using EEGLAB ----
     %disp('Using imported events...');
    %events = readtable(fullfile(data_dir, this_events));
    %events.latency = events.latency + 1; % add 1 as EMSE starts at 0
    %EEG.event = table2struct(events); % inserts outside events
    %else
     % File does not exist...load standard file
    % Loads in raw data using EEGLAB ----
    %disp('Using embedded events...');
    %end
    
    % Remove externals that are not being used ----
    EEG = pop_select(...
        EEG,...
        'rmchannel',{'EXG1', 'EXG2', 'EXG3','EXG4','EXG5','EXG6','EXG7','EXG8'}...
        );
    
    % Checks to see if there are more than 64 channels in the recording----
    % first creates a vector of the chan names   
    N = 32;
    A_chans = cell(1, N);
    B_chans = A_chans;
    for j=1:N
        A_chans{j} = strcat('A', num2str(j));
        B_chans{j} = strcat('B', num2str(j));
    end
    chans_to_keep = [A_chans B_chans]; % here is the vector of channel names
    
    % keeps these chans only
    if EEG.nbchan > 64
        EEG = pop_select(EEG, 'channel', chans_to_keep);
    else
        disp('only 64 channels detected; not removing any...');
    end
        
    % Configuring channel locations ----
    % loads in ELP
    %eloc = readlocs( this_elp ); % reads in elp chan locations
    eloc = readlocs(this_elp, ...
        'filetype', 'custom', ... % uses custom file I created
        'skiplines', 1, ... % skips header
        'format', {'channum','labels','X','Y','Z'});
    EEG.chanlocs = eloc; % adds chan locs
    
    % sets A1 as ref because it was chosen upon import
    EEG = pop_chanedit(EEG, 'setref', {'1:64' 'Fp1'}); 
    
    % Downsamples to 256Hz ----
    EEG = pop_resample(EEG, 256);
    
    % Re-references data to 
    %EEG = pop_reref( EEG, [24 61] ,'keepref','on'); % linked mastoids
    EEG = pop_reref( EEG, []); % common average reference
    
    % Removing DC offset by subtracting the mean signal from each electrode
    EEG = pop_rmbase(EEG, [], []);
    
    % Highpass filter at 1Hz 
    % -6dB @ 1Hz, 425 point highpass, 2Hz transition band width
    EEG = pop_eegfiltnew(EEG, 'locutoff', 2, 'plotfreqz', 0);

    % Cleanline ----
    % Removing electrical line noise @ 60 Hz
    EEG = pop_cleanline(EEG, 'bandwidth', 2, 'chanlist', [1:EEG.nbchan],...
        'computepower', 1, 'linefreqs', 60, 'normSpectrum', 0, ...
        'p', 0.01, 'pad', 2, 'plotfigures' , 0, 'scanforlines', 1, ...
        'sigtype', 'Channels', 'tau', 100, 'verb', 1, 'winsize', ...
        4,'winstep',1);
    
    % Renames dataset ----
    dataset_name = strcat(this_ss, '-prepro');
    EEG = pop_editset(EEG, 'setname', dataset_name, 'run', []);

    % Saves out preprocessed data for inspection ----
    EEG = pop_saveset(EEG, 'filename', dataset_name, 'filepath', output_dir);
    
    eeglab redraw % redraws to GUI for convenience

end