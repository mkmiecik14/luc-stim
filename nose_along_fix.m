% nose along +Y
% Matt Kmiecik
% 04 May 2024

% Purpose: this script details the work done to save out the channel
% coordinates that are rotated from having the nose along +X to +Y. This
% was done because I initially used the coordiantes provided from Biosemi's
% website that habe the nose along +Y, but EEGLAB assumes they are along +X
% This has unintended consequences for ICA because IClabel uses topo maps
% to guide labelling, therefore causing it to underperform because the 
% topomaps were incorrect. Because nose along cannot be changed 
% programatically, a new chan locations file will be imported and 
% substituted before IC labelling

workspace_prep

this_ss = '1461831842022_2'; % subject number does not really matter
this_ss_path = dir(fullfile(output_dir, strcat(this_ss, '-ica.set')));
this_ss_name = this_ss_path.name;
        
% Loads in data using EEGLAB ----
EEG = pop_loadset('filename',this_ss_name,'filepath', this_ss_path.folder);
eeglab redraw

% from here I went in manually into Edit > Channel locations > and changed
% the nose along +X to +Y. I saved out a *.ced file as well

% saves into variables
chan_info = EEG.chaninfo;
chan_locs = EEG.chanlocs;

% saves out as matlab struct ( just in case they are needed )
save("doc\chan_info_nose_along_fixed.mat", 'chan_info');
save("doc\chan_locs_nose_along_fixed.mat", 'chan_locs'); 
