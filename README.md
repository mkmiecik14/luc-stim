# luc-stim

Code for processing and analyzing resting state EEG data from the CANlab

The main directory for the project looks like this:

```
luc-stim
├───data
├───doc
├───matlab-toolboxes
├───output
└───src
```

The contents of this github repo are of `src/`. Therefore, paths are relative to `luc-stim/`

`data/` contains all raw data (e.g., BDFs, ELPs, etc.)
`doc/` contains documentation and notes. Importantly, the Excel worksheet `ss-info.xlsx` is crucial to the pipeline
`matlab-toolboxes/` contains all the MATLAB toolboxes for the project. These must be added to the path within MATLAB*
`output/` collects all output produced from scripts
`src/` contains all scripts and functions for the project and is the github repo

* can be located elsewhere, like in `C:\Documents`, but it is good practice to keep these here for version control.

# EEG Pipeline

The script `workspace_prep.m` is run at the beginning of each processing step and prepares the MATLAB work space. **The paths in this script must be updated to the machine.**

The "batch" sheet in `docs/ss-info.xlsx` must be updated with the desired participant(s).  

## Step 1 - Preprocessing

This step preprocesses the EEG data by importing BDFs into EEGLAB, removing unnecessary channels, configuring channel locations, downsampling, computing common average reference, removing DC offset, high pass filtering, removing 60Hz line noise, and saves out dataset into `output/` for visual inspection.

Visually inspect data after this step and note any bad channels. Bad channels must be entered into the second column in the "batch worksheet within `doc/ss-info.xlsx` inside of square brackets according to index. For example, suppose that channel Fp1 and Pz were deemed bad and in need of interpolation. The `bad_chans` column would be [1 31]

## Step 2- ICA

ICA is run on participants using the `eeg_ica.m` script in EEGLAB. Bad channels are interpolated and the rank of the matrix is adjusted accordingly prior to ICA decomposition.

## Step 3 - FFT

Running the script `eeg_spect.m` will perform artifact correction by removing ICs. Then it will perform individual alpha frequency calculations for each block, both peak alpha frequency and center of gravity estimates as well as whole scalp estimates. Next, pink and white noise estimation occurs and is subtracted from the signal (currently saved separately).

Plots of block-wise broadband frequency and pink/white noise estimation are saved in `output/` for blocks that are available. All results are saved in a matlab structure named `ssnumber-spec-res.mat` for import into R or your favorite analysis software.

Note: for pink and white noise estimation to be reflected in IAF calculations, this script will need to be modified.

## Step 4 - Preprocess MATLAB files to R

The next step is to take the spectral results MATLAB structure file `ssnumber-spec-res.mat` and get these files ready for statistical processing in R. This is accomplished by `prepro-to-r.R`.
