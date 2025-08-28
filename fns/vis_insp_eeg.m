function EEG = vis_insp_eeg(ss)

% Function for visualizing EEG signal post ICA
% Matt Kmiecik
% Started 06 APRIL 2024
% USAGE EXAMPLE: EEG = viz_eeg_post_ica('1461831842001_1'); eeglab redraw

workspace_prep % Prepares workspace
    
% Creating variables ----
this_ss_path = dir(fullfile(output_dir, strcat(ss, '-prepro.set')));
this_ss_name = this_ss_path.name;
    
% Loads in data using EEGLAB ----
EEG = pop_loadset('filename',this_ss_name,'filepath', this_ss_path.folder);

eeglab redraw % redraws

% plots
eegplot(EEG.data, 'spacing', 100, 'winlength', 30, 'events', EEG.event);

end