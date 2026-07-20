clear;
close all;
clc;

rootDir = fileparts(mfilename('fullpath'));
projectDir = fileparts(rootDir);
addpath(rootDir);
addpath(projectDir, '-end');

whale_file_local = fullfile(rootDir, 'whale_7350_1.mat');
whale_file_parent = fullfile(projectDir, 'whale_7350_1.mat');
if exist(whale_file_local, 'file')
    load(whale_file_local);
elseif exist(whale_file_parent, 'file')
    load(whale_file_parent);
else
    error('whale_7350_1.mat was not found in %s or %s.', rootDir, projectDir);
end

fs = 7350;
signal = double(data(:));
signal = filter([1 -1], 1, signal);

win_len = 581;
time_points = 380;
fre_len = 824;
chirp_len = 331;
threshold = 0.029;

pcct_params = struct();
pcct_params.use_local_predict_score = true;
pcct_params.use_projection_amp_comp = true;

[tfr_final, tt, freq] = PCCT( ...
    signal, fs, win_len, time_points, fre_len, chirp_len, ...
    threshold, pcct_params);

tfr_plot = abs(tfr_final);
peak_value = max(tfr_plot(:));
if peak_value <= 0
    peak_value = eps;
end
tfr_plot = tfr_plot / peak_value;
display_floor_db = -60;
tfr_plot_db = 20 * log10(tfr_plot.' + eps);
tfr_plot = (max(tfr_plot_db, display_floor_db) - display_floor_db) / ...
    abs(display_floor_db);

signal_duration = length(signal) / fs;
pcct_plot_tfr(tt, freq, tfr_plot, [0 1], ...
    [0 signal_duration], [0 3500], ...
    'Time [s]', 'Frequency [Hz]', 0.163);
