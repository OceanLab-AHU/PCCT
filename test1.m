clear;
close all;
clc;

rootDir = fileparts(mfilename('fullpath'));
projectDir = fileparts(rootDir);
addpath(rootDir);
addpath(projectDir, '-end');

fs = 100;
t_end = 10;
t = (0:1 / fs:t_end - 1 / fs)';

if1 = 25 + 15 * sin(0.311 * pi * (1.1 * t + 1.9));
if2 = 25 - 15 * sin(0.311 * pi * (1.1 * t + 1.9));

phi1 = 2 * pi * cumtrapz(t, if1);
phi2 = 2 * pi * cumtrapz(t, if2);

signal_clean = exp(1j * phi1) + exp(1j * phi2);
snr_db = 0;
if exist('awgn', 'file') == 2
    signal = awgn(signal_clean, snr_db, 'measured');
else
    signal = add_awgn_local(signal_clean, snr_db);
end

threshold = 0.4;
win_len = 121;
signal_len = length(signal);
fre_len = signal_len;
time_points = signal_len;
chirp_len = 111;

pcct_params = struct();
pcct_params.use_local_predict_score = true;
pcct_params.use_projection_amp_comp = true;

[tfr_final, tt, freq] = PCCT( ...
    signal, fs, win_len, time_points, fre_len, chirp_len, ...
    threshold, pcct_params);

tfr_plot = abs(tfr_final).';
peak_value = max(tfr_plot(:));
if peak_value <= 0
    peak_value = eps;
end

pcct_plot_tfr(tt, freq, tfr_plot, [0 peak_value], ...
    [0 t_end], [0 fs / 2], ...
    'Time [s]', 'Frequency [Hz]', 0.183);

function y = add_awgn_local(x, snr_db)
x = x(:);
sig_power = mean(abs(x).^2);
noise_power = sig_power / (10^(snr_db / 10));
if isreal(x)
    noise = sqrt(noise_power) * randn(size(x));
else
    noise = sqrt(noise_power / 2) * (randn(size(x)) + 1j * randn(size(x)));
end
y = x + noise;
end
