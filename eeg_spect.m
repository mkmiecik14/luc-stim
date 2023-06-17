% EEG ICA cleaning and Spectral Decomposition Step 3
% Matt Kmiecik
% Started 03 JUNE 2023

workspace_prep % Prepares workspace

% Preallocation ----
ss = string({RAW{2:size(RAW,1),1}}); % ss strings
num_iters = length(ss);         % number of participants in this batch
i=1;                            % for testing purposes
csd_switch = 0;                 % 1 == CSD will be computed
plot_switch = 1;                % 1 == PSD plots will be saved
noise_plot_switch = 1;          % 1 == Pink&White Noise plots saved

% For spectral decomposition settings
wsize = 4;      % FFT window size in seconds
olap = 2;       % FFT window overlap in seconds

% For peak alpha frequency (PAF) and center of gravity (COG) calculations
% see: https://github.com/corcorana/restingIAF/blob/master/tutorial/tutorial.m
cmin = 3;           % minimum number of channel estimates required for 
                    % cross-channel averages
fRange = [1 40];    % spectral range (set to filter passband)
w = [7.5 13];       % alpha peak search window (Hz); from McClain et al. 2022
Fw = 23;            % SGF frame width (Fw * freq resoltion = window)
poly = 5;           % SGF polynomial order (same as described in Corcoran et al., 2017)

for i = 1:num_iters
    
    % Creating variables ----
    this_ss = ss(i);
    this_ss_path = dir(fullfile(output_dir, strcat(this_ss, '-ica.set')));
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
    
    eeglab redraw    
   
    if csd_switch == 1
    
        % Compute surface laplacian spatial filter ----
        % Calculating Surface Laplacian via CSD toolbox functions
        % http://psychophysiology.cpmc.columbia.edu/software/CSDtoolbox/tutorial.html#PrepareInput
        chan_mont = cell(EEG.nbchan,1); % initializes cell array
        % Fills cell array with electrode labels
        for j = 1:length(chan_mont)
            chan_mont(j) = cellstr(EEG.chanlocs(j).labels);
        end
        % Derives spherical coordinates via CSD toolbox fx ExtractMontage()
        csd_mont = ExtractMontage('10-5-System_Mastoids_EGI129.csd', chan_mont);
        % To view: MapMontage(csd_mont)
        [G, H] = GetGH(csd_mont); % Calculates G and H matrices

        % Applies surface laplacian to EEG data
        % lambda left at default, 10 = cm head size, so units are microvolt/cm^2
        EEG.data = CSD(EEG.data, G, H, 1.0e-5, 10); 
    
    else
        disp('CSD skipped....');
    end
    
    % Spectral decomposition ----
    
    % Epoching
    % selecting the stimulation blocks (60 second epochs)
    blocks = {'101','103','105','107', '109', '203', '205', '207', '209'};
    %blocks_end = {'S211' 'S202' 'S203' 'S214' 'S205' 'S216' 'S217' 'S208'};
    
    % preallocates arrays
    N = EEG.srate*wsize; % number of data points in FFT window
    poss_freqs = N/2 + 1; % possible frequencies (freq bins = EEG.srate/N)
    this_spectra = NaN(EEG.nbchan, poss_freqs, length(blocks));
    this_spectra_psd = this_spectra; % for conversion to psd (no decibels)
    this_spectra_amp = this_spectra; % for conversion to amplitude
    this_freqs = NaN(poss_freqs, 1, length(blocks));
    this_paf = NaN(EEG.nbchan, length(blocks)); % peak alpha freq
    this_cog = NaN(EEG.nbchan, length(blocks)); % center of gravity
    
    % grand average PAF and COG (rows are blocks, cols are PAF and COG)
    this_iaf = NaN(length(blocks), 2); 
    
    % In order for this script to accomodate rejected epochs, the
    % epoch start stop must be determined by a 
    % starting trigger (e.g., S111) and its ending e.g., (S211)
    
    % Getting start and stop latencies of the stimulation blocks
%     start_times = NaN(1, length(blocks)); % initializes vector
%     end_times = start_times; % initializes vector
%     for k = 1:length(blocks)
%         % Retrieving the latency (in seconds) of each block start
%         start_times(1,k) = (eeg_getepochevent(EEG, blocks(k), [], 'latency'))/1000;
%         % Retrieving the latency (in seconds) of each block end
%         end_times(1,k) = (eeg_getepochevent(EEG, blocks_end(k), [], 'latency'))/1000;
%     end
    
    % Calculates the duration of each block (in order)
    %block_durations = round(end_times - start_times);
    
    j=1; % for testing purposes
    
    for j = 1:length(blocks)
        
        this_epoch = [0 60]; % Sets the epoch duration
        
        try % This will run if the block exists
            % Selects blocks (in order)
            this_EEG = pop_epoch(EEG,blocks(j),this_epoch,'epochinfo','yes');
            % Spectral decomposition here
            [this_spectra(:,:,j), this_freqs(:,:,j)] = spectopo(...
                this_EEG.data(:,:), ... 
                0, ... % frames per epoch (default 0 = data length)
                this_EEG.srate, ... % sampling rate
                'winsize', wsize*this_EEG.srate, ... % window size
                'overlap', olap*this_EEG.srate, ... % overlap size
                'plot','off'... % toggles plot
                );
            
            % obtains frequency resolution
            freq_res = this_freqs(2,1,j)-this_freqs(1,1,j);
            
            % Converts from dB to PSD (microvolts^2/Hz)
            this_spectra_psd(:,:,j) = 10.^(this_spectra(:,:,j)/10);
            
            % Converts from PSD (microvolts^2/Hz) to amplitdue (microvolts^2)
            % I decided to keep things in terms of PSD (microvolts^2/Hz)
            %this_spectra_amp(:,:,j) = freq_res.*this_spectra_psd(:,:,j);
             
            % Peak alpha frequency (PAF) and center of gravity (COG)
            [pSum, pChans, f] = restingIAF(...
                this_EEG.data,... % EEG data
                this_EEG.nbchan,... % number of channels
                cmin,... % min number of channels to grand average PAF
                fRange,... % spectral range
                this_EEG.srate,... % sampling rate
                w,... % alpha peak search window
                Fw,... % window size
                poly,... % polynomial order
                'taper', 'hamming',... % type of taper (default)
                'tlen', N,...% so that window size remains the same as spectral decomp (50% overlap by default)
                'mdiff', .2); % peak must be 20% larger than surrounding (default)  
            
            % fills in preallocated matrices
            this_paf(:,j) = [pChans.peaks]; % peak alpha freq per chan
            this_cog(:,j) = [pChans.gravs]; % center of gravity per chan
            this_iaf(j,:) = [pSum.paf pSum.cog]; %IAF (weighted by Q; see Corcoran et al., 2017)
            
            % Plots for troubleshooting (if needed)
            if plot_switch == 1
                figure("Visible","off");
                pop_spectopo(this_EEG,1,[],'EEG','freq',[4 8 12 25 30],'freqrange',[0 75],'electrodes','on');
                saveas(gcf,...
                    fullfile(output_dir,...
                    strcat(num2str(this_ss),...
                    '-',blocks{j},'.png')));
                close; % closes figure
            else
                % plots not saved
            end
        
        catch % if the block is missing, then the matrix is filled with NaN
            this_spectra(:,:,j) = NaN; % fills with missing values
            this_freqs(:,:,j) = NaN; % fills with missing values
     
        end
                           
    end
    
    % Calculating pink and white noise via PaWNextra.m
    % See Barry & De Blasio (2021)
    %
    % ! The minimum frequency to begin PaWNextra is important !
    % given that these data were highpass filtered (see
    % prepro_resting_state.m) we should start estimating pink&white noise
    % at the edge of the filter (i.e., if data were highpass filtered 
    %  -6dB @ 1Hz, 425 point highpass, 2Hz transition band width, then the
    %  minimum frequency to start noise estimation would be 2 Hz
    min_freq = 2; % minimum frequency to begin noise estimation
    % This loop was developed to determine the frequency vector in the case
    % of a missing/rejected first block(s)
    for m=1:size(this_freqs,3)
        if isnan(this_freqs(:,1,m))
            % do nothing as this block's frequency vector is nan
        else
            non_dc_index = find(this_freqs(:,:,m) >= min_freq);
            freq_vector = this_freqs(non_dc_index,:,m);
            break % no need to evaluate further as frequency vector is determined
        end
    end

    this_spectra_psd_t = permute(this_spectra_psd(:,non_dc_index,:), [2 1 3]); % transposes
    [iterations, PN, PN_slope, WN, n_sols, adjust_value, error_value, FL_chan] = ...
            PaWNextra(this_spectra_psd_t,... % broadband spectra freqs x chans x participants)
            freq_vector,... % column vector of freq bins
            -.1,... % error threshold (E=-.001 is the default)
            1000); % iterations (I believe this is the default)
    
    % Subtracting pink and white noise from broadband spectra
    this_corrected = NaN(size(this_spectra_psd_t)); % initialization
    for slice=1:size(this_spectra_psd_t,3)
        this_WN = repmat(WN(slice,:), [length(freq_vector),1]); % preps
        % Pink and White Noise correction HERE!
        this_corrected(:,:,slice) = this_spectra_psd_t(:,:,slice) - PN(:,:,slice) - this_WN;
        % PLOTS RESULTS AND SAVES (if requested):
        if noise_plot_switch == 1
            if ~isnan(this_spectra_psd_t(1,1,slice))
                figure("Visible", "off"); % plots the inital, pink, and observed spectra
                axes('FontSize', 6)
                t=tiledlayout(8,8,'TileSpacing','compact');
                for m=1:EEG.nbchan % was 30
                    nexttile
                    % Initial spectra
                    plot(freq_vector, this_spectra_psd_t(:,m,slice), 'Color', [.7 .7 .7])
                    hold on
                    % Pink Noise
                    plot(freq_vector, PN(:,m,slice), '-m')
                    hold on
                    % Observed spectra
                    plot(freq_vector, this_corrected(:,m,slice)', '-k')
                    axis([0 15 0 Inf])
                    title(EEG.chanlocs(m).labels,'FontSize',6)
                end
                title(t, strcat("Portcode", blocks{slice}))
                xlabel(t, 'Hz')
                ylabel(t, 'PSD [uV^2/Hz]') % after surface Laplacian, this is (uV/cm^2)^2/Hz (sort of like PSD)
                plot_name = strcat(num2str(this_ss),...
                        '-',blocks{slice},'-pink.pdf');
                set(gcf,'units','inches','Position',[2, 2, 8, 6])
                exportgraphics(gcf,fullfile(output_dir, plot_name),'ContentType','vector');
                close; % closes figure down
            else
                % figure is not drawn/saved because block is NaN
            end
        else
            % figures not drawn/saved because plots were not requested
        end
    end
    
    % Saving out results ----
    % combines into one variable
    % stimulation results are stored with each index being the
    % blocks in the following order:
    % 101 - eyes open 182_IAF
    % 103 - eyes open 183_PRE_RS
    % 105 - eyes closed 183_PRE_RS
    % 107 - eyes open 183_PRE_RS
    % 109 - eyes closed 183_PRE_RS
    % 203 - eyes open 184_POST_RS
    % 205 - eyes closed 184_POST_RS
    % 207 - eyes open 184_POST_RS
    % 209 - eyes closed 184_POST_RS
    spec_res.spectra   = this_spectra; % 3D mat of spectra (in dB)
    spec_res.freqs     = this_freqs;   % 3D mat of freqs bins
    spec_res.paf       = this_paf;     % 2D mat of PAF per chan per block
    spec_res.cog       = this_cog;     % 2D mat of COG per chan per block
    spec_res.iaf       = this_iaf;     % 2D mat of IAF per block
    spec_res.psd       = this_spectra_psd; % converted this_spectra from dB to PSD (uV^2/Hz)
    spec_res.corrected = this_corrected; % Pink and white noise are subtracted from PSD
    spec_res.freqvec   = freq_vector; % frequency vector used for Pink&White noise calc

    % Saving out all data ----
    spec_outname = strcat(this_ss, '-spec-res.mat');
    save(fullfile(output_dir, spec_outname), 'spec_res'); % saves out as matlab struct
end