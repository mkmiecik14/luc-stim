function EEG = viz_eeg_post_ica(ss)

% Function for visualizing EEG signal post ICA
% Matt Kmiecik
% Started 06 APRIL 2024
% USAGE EXAMPLE: EEG = viz_eeg_post_ica('1461831842001_1'); eeglab redraw

workspace_prep % Prepares workspace
    
% Creating variables ----
this_ss_path = dir(fullfile(output_dir, strcat(ss, '-ica.set')));
this_ss_name = this_ss_path.name;
    
% Loads in data using EEGLAB ----
EEG = pop_loadset('filename',this_ss_name,'filepath', this_ss_path.folder);

% Labels ICs for rejection ----
EEG = pop_iclabel(EEG, 'default');
EEG = pop_icflag(EEG, ...
    [NaN NaN;...    % brain
    0.8 1;...       % muscle (> 80% probability will reject components)
    0.8 1;...       % eye (> 80% probability will reject components)
    NaN NaN;...     % heart
    NaN NaN;...     % line noise
    NaN NaN;...     % channel noise
    NaN NaN...      % other
    ]);

% Removes artifactual ICs
this_reject = find(EEG.reject.gcompreject);
EEG = pop_subcomp(EEG, this_reject, 0);

end