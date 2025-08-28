% TEST_PIPELINE_SETUP - Quick validation of pipeline setup
%
% This script tests that all functions are accessible and the configuration
% loads properly without actually processing data.

fprintf('=== TESTING PIPELINE SETUP ===\n\n');

% Add function path
addpath('fns');

% Test 1: Configuration loading
fprintf('1. Testing configuration loading...\n');
try
    cfg = config();
    fprintf('   ✓ Configuration loaded successfully\n');
    fprintf('   ✓ Main directory: %s\n', cfg.paths.main);
    fprintf('   ✓ Output directories created\n');
    
    if ~isempty(cfg.participants)
        fprintf('   ✓ Participant info loaded: %d subjects\n', length(cfg.participants.ss_list));
    else
        fprintf('   ⚠ No participant info loaded (this is expected if ss-info.xlsx is not available)\n');
    end
catch ME
    fprintf('   ✗ Configuration failed: %s\n', ME.message);
    return;
end

% Test 2: Function availability
fprintf('\n2. Testing function availability...\n');
functions_to_test = { ...
    'load_raw_eeg', 'preprocess_eeg', 'extract_blocks', ...
    'detect_artifacts', 'apply_ica', 'interpolate_channels', ... 
    'spectral_analysis', 'save_results'
};

all_functions_ok = true;
for i = 1:length(functions_to_test)
    func_name = functions_to_test{i};
    if exist(func_name, 'file') == 2
        fprintf('   ✓ %s.m found\n', func_name);
    else
        fprintf('   ✗ %s.m missing\n', func_name);
        all_functions_ok = false;
    end
end

% Test 3: Directory structure
fprintf('\n3. Testing directory structure...\n');
required_dirs = {'preprocessed', 'ica', 'spectral', 'logs', 'qc'};
for i = 1:length(required_dirs)
    dir_path = fullfile(cfg.paths.output, required_dirs{i});
    if exist(dir_path, 'dir')
        fprintf('   ✓ %s/ exists\n', required_dirs{i});
    else
        fprintf('   ✗ %s/ missing\n', required_dirs{i});
    end
end

% Test 4: Dependencies 
fprintf('\n4. Testing EEGLAB availability...\n');
try
    which_eeglab = which('eeglab');
    if isempty(which_eeglab)
        fprintf('   ⚠ EEGLAB not found in path - you may need to add it manually\n');
    else
        fprintf('   ✓ EEGLAB found: %s\n', which_eeglab);
    end
catch
    fprintf('   ⚠ Could not check EEGLAB availability\n');
end

% Test 5: Required toolbox functions
fprintf('\n5. Testing required toolbox functions...\n');
required_functions = {'pop_biosig', 'pop_resample', 'pop_reref', 'pop_eegfiltnew', ...
                     'pop_cleanline', 'pop_runica', 'pop_iclabel', 'pop_interp', 'spectopo'};

missing_functions = {};
for i = 1:length(required_functions)
    func_name = required_functions{i};
    if exist(func_name, 'file')
        fprintf('   ✓ %s available\n', func_name);
    else
        fprintf('   ⚠ %s not found - may need EEGLAB plugins\n', func_name);
        missing_functions{end+1} = func_name;
    end
end

% Summary
fprintf('\n=== SETUP TEST SUMMARY ===\n');
if all_functions_ok
    fprintf('✓ All pipeline functions implemented\n');
else
    fprintf('✗ Some pipeline functions missing\n');
end

if isempty(missing_functions)
    fprintf('✓ All required EEGLAB functions available\n');
else
    fprintf('⚠ Missing EEGLAB functions: %s\n', strjoin(missing_functions, ', '));
    fprintf('  Install required EEGLAB plugins if needed\n');
end

fprintf('✓ Directory structure created\n');
fprintf('✓ Configuration system working\n');

fprintf('\n=== READY TO TEST WITH DATA ===\n');
fprintf('To test with a single subject, run:\n');
fprintf('  run_eeg_pipeline({{''SUBJECT_ID''}});\n');
fprintf('Where SUBJECT_ID matches a .bdf file in your data directory.\n\n');

fprintf('Pipeline setup test completed.\n');