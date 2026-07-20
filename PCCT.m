function [TFR, time_grid, freq_grid] = PCCT( ...
    signal, fs, win_len, time_points, fre_len, chirp_len, threshold, ...
    pcct_params, soft_mask_factor)
if nargin < 8 || isempty(pcct_params)
    pcct_params = default_pcct_params();
end
if nargin < 9 || isempty(soft_mask_factor)
    soft_mask_factor = 1;
end

pcct_params = fill_pcct_defaults(pcct_params);

[TFR, time_grid, freq_grid] = adaptive_decay_main( ...
    signal, fs, win_len, time_points, fre_len, chirp_len, threshold, ...
    pcct_params, soft_mask_factor);
end

function [TFR, time_grid, freq_grid] = adaptive_decay_main( ...
    signal, fs, win_len, time_points, fre_len, chirp_len, threshold, ...
    pcct_params, soft_mask_factor)

signal = signal(:).';
signal_len = length(signal);
fre_half = round(fre_len / 2);

chirp_len = chirp_len + 1 - rem(chirp_len, 2);
win_len = win_len + 1 - rem(win_len, 2);
win_half = (win_len - 1) / 2;

chirp_ext = chirp_len * 2 - 1;
c_end = fs * fs / win_len;

freq_grid_ext = (0:(fre_len - 1)) * fs / fre_len;
freq_grid = freq_grid_ext(1:fre_half);
time_grid = (1 / time_points:1 / time_points:1) * signal_len / fs;
chirp_grid = (-1:2 / (chirp_len - 1):1) * c_end;
chirp_grid_ext = (-2:2 / (chirp_len - 1):2) * c_end;

sig_tfcr = zeros(time_points, fre_len, chirp_len);

alpha = 3;
if exist('gausswin', 'file') == 2
    win = gausswin(win_len, alpha);
else
    n = (-win_half:win_half).';
    sigma_sample = win_len / 6;
    win = exp(-(n.^2) / (2 * sigma_sample^2));
end
win = win / sum(win);

tau_index = round(time_grid * fs);

if isreal(signal)
    signal = analytic_signal_fft(signal);
end

signal_ext = [zeros(1, win_half), signal, zeros(1, win_half)];
win_time = (-win_half:win_half) / fs;

sig_mat = signal_ext( ...
    repmat((0:win_len - 1), [time_points, 1, chirp_len]) + ...
    repmat(tau_index.', [1, win_len, chirp_len]));

win_mat = repmat(win.', [time_points, 1, chirp_len]);

kernel = zeros(1, win_len, chirp_len);
kernel(1, :, :) = exp(-1j * pi * (win_time).^2.' * chirp_grid);
kernel_mat = repmat(kernel, [time_points, 1, 1]);

sig_tfcr_a = sig_mat .* win_mat .* kernel_mat;
sig_tfcr(:, 1:win_len, :) = sig_tfcr_a;

sig_tfcr_c = circshift(sig_tfcr, -win_half, 2);
sig_tfcr_c = fft(sig_tfcr_c, [], 2);
sig_tfcr = sig_tfcr_c(:, 1:fre_half, :);
sig_tfcr_origin = sig_tfcr;

ref_tfcr = zeros(fre_len, chirp_ext);
ref_win_mat = repmat(win, [1, chirp_ext]);
ref_kernel_mat = exp(-1j * pi * (win_time).^2.' * chirp_grid_ext);

ref_tfcr_a = ref_win_mat .* ref_kernel_mat;
ref_tfcr(1:win_len, :) = ref_tfcr_a;

ref_tfcr_c = circshift(ref_tfcr, -win_half, 1);
ref_tfcr_c = fft(ref_tfcr_c, [], 1);
lct_ref = fftshift(ref_tfcr_c, 1);

freq_center = round((fre_len + 1) / 2) * ones(1, time_points);
chirp_center = (chirp_ext + 1) / 2 * ones(1, time_points);

mf = 0:fre_half - 1;
mc = zeros(1, 1, chirp_len);
mc(1, 1, :) = (0:chirp_len - 1) * fre_len;

mfc = repmat(mf, time_points, 1, chirp_len) + ...
    repmat(mc, time_points, fre_half, 1);

mid_c = (chirp_len + 1) / 2;
center_slice = sig_tfcr(:, :, mid_c);
signal_energy = sum(abs(center_slice(:)).^2) + eps;

ratio = 1;
iter_count = 0;

init_points = pcct_params.init_points;
freq_range = pcct_params.freq_range;
chirp_range = pcct_params.chirp_range;
use_local_predict_score = pcct_params.use_local_predict_score;
max_freq_ridges = {};
max_chirp_ridges = {};
max_chirp_ridges_observed = {};

while ratio > threshold
    iter_count = iter_count + 1;

    abs_sig_tfcr = abs(sig_tfcr);
    global_tfcr_peak = max(abs_sig_tfcr(:)) + eps;

    best_energy = -Inf;
    best_freq_ridge = zeros(1, time_points);
    best_chirp_ridge = zeros(1, time_points);
    best_chirp_ridge_observed = zeros(1, time_points);
    best_mode_ridge = zeros(1, time_points);

    search_points = round(linspace( ...
        time_points / (init_points + 1), ...
        time_points - time_points / (init_points + 1), ...
        init_points));
    search_points = max(2, min(time_points - 1, search_points));
    search_points = unique(search_points);

    for init_idx = 1:numel(search_points)
        k0 = search_points(init_idx);

        slice0 = squeeze(abs_sig_tfcr(k0, :, :));
        [~, idx0] = max(slice0(:));
        [f0, c0] = ind2sub([fre_half, chirp_len], idx0);

        temp_freq_ridge = zeros(1, time_points);
        temp_chirp_ridge = zeros(1, time_points);
        temp_chirp_ridge_observed = zeros(1, time_points);
        temp_mode_ridge = zeros(1, time_points);

        temp_freq_ridge(k0) = f0;
        temp_chirp_ridge(k0) = c0;
        temp_chirp_ridge_observed(k0) = c0;
        temp_idx = sub2ind([time_points, fre_half, chirp_len], k0, f0, c0);
        temp_mode_ridge(k0) = sig_tfcr(temp_idx);

        for nt = k0 + 1:time_points
            prev_f = temp_freq_ridge(nt - 1);
            prev_c = temp_chirp_ridge(nt - 1);

            f_set = max(1, prev_f - freq_range):min(fre_half, prev_f + freq_range);
            c_set = max(1, prev_c - chirp_range):min(chirp_len, prev_c + chirp_range);

            patch = reshape(abs_sig_tfcr(nt, f_set, c_set), ...
                numel(f_set), numel(c_set));
            [~, observed_idx] = max(patch(:));
            [~, observed_c_off] = ind2sub(size(patch), observed_idx);
            if use_local_predict_score
                dt_local = time_grid(nt) - time_grid(nt - 1);
                score = local_predict_score(patch, f_set, c_set, prev_f, prev_c, ...
                    dt_local, freq_grid, chirp_grid, freq_range, chirp_range, ...
                    global_tfcr_peak);
            else
                score = patch;
            end
            [~, local_idx] = max(score(:));
            [f_off, c_off] = ind2sub(size(score), local_idx);

            new_f = f_set(f_off);
            new_c = c_set(c_off);

            temp_freq_ridge(nt) = new_f;
            temp_chirp_ridge(nt) = new_c;
            temp_chirp_ridge_observed(nt) = c_set(observed_c_off);

            temp_idx = sub2ind([time_points, fre_half, chirp_len], nt, new_f, new_c);
            temp_mode_ridge(nt) = sig_tfcr(temp_idx);
        end

        for nt = k0 - 1:-1:1
            next_f = temp_freq_ridge(nt + 1);
            next_c = temp_chirp_ridge(nt + 1);

            f_set = max(1, next_f - freq_range):min(fre_half, next_f + freq_range);
            c_set = max(1, next_c - chirp_range):min(chirp_len, next_c + chirp_range);

            patch = reshape(abs_sig_tfcr(nt, f_set, c_set), ...
                numel(f_set), numel(c_set));
            [~, observed_idx] = max(patch(:));
            [~, observed_c_off] = ind2sub(size(patch), observed_idx);
            if use_local_predict_score
                dt_local = -(time_grid(nt + 1) - time_grid(nt));
                score = local_predict_score(patch, f_set, c_set, next_f, next_c, ...
                    dt_local, freq_grid, chirp_grid, freq_range, chirp_range, ...
                    global_tfcr_peak);
            else
                score = patch;
            end
            [~, local_idx] = max(score(:));
            [f_off, c_off] = ind2sub(size(score), local_idx);

            new_f = f_set(f_off);
            new_c = c_set(c_off);

            temp_freq_ridge(nt) = new_f;
            temp_chirp_ridge(nt) = new_c;
            temp_chirp_ridge_observed(nt) = c_set(observed_c_off);

            temp_idx = sub2ind([time_points, fre_half, chirp_len], nt, new_f, new_c);
            temp_mode_ridge(nt) = sig_tfcr(temp_idx);
        end

        temp_energy = sum(abs(temp_mode_ridge).^2);

        if temp_energy > best_energy
            best_energy = temp_energy;
            best_freq_ridge = temp_freq_ridge;
            best_chirp_ridge = temp_chirp_ridge;
            best_chirp_ridge_observed = temp_chirp_ridge_observed;
            best_mode_ridge = temp_mode_ridge;
        end
    end

    best_freq_ridge = max(1, min(fre_half, round(best_freq_ridge)));
    best_chirp_ridge = max(1, min(chirp_len, round(best_chirp_ridge)));
    best_chirp_ridge_observed = max(1, min(chirp_len, round(best_chirp_ridge_observed)));
    smooth_win = 7;
    best_freq_ridge_s = smooth_ridge_indices(best_freq_ridge, smooth_win, fre_half);
    best_chirp_ridge_s = smooth_ridge_indices(best_chirp_ridge, smooth_win, chirp_len);

    max_freq_ridges{end + 1} = best_freq_ridge_s; %#ok<AGROW>
    max_chirp_ridges{end + 1} = best_chirp_ridge_s; %#ok<AGROW>
    max_chirp_ridges_observed{end + 1} = best_chirp_ridge_observed; %#ok<AGROW>

    for nt = 1:time_points
        f_idx = best_freq_ridge(nt);
        c_idx = best_chirp_ridge(nt);

        fdist = freq_center(nt) - f_idx + 1;
        cdist = chirp_center(nt) - c_idx + 1;

        mtau = fdist + (cdist - 1) * fre_len;
        m1b = repmat(mtau, 1, fre_half, chirp_len) + mfc(nt, :, :);
        m1b = max(1, min(numel(lct_ref), m1b));

        sig_tfcr_k = repmat(best_mode_ridge(nt), 1, fre_half, chirp_len) .* lct_ref(m1b);
        sig_tfcr(nt, :, :) = sig_tfcr(nt, :, :) - soft_mask_factor * sig_tfcr_k;
    end

    center_slice_res = sig_tfcr(:, :, mid_c);
    residual_energy = sum(abs(center_slice_res(:)).^2);
    ratio = residual_energy / signal_energy;
end

TFR = render_tfr_from_ridges( ...
    abs(sig_tfcr_origin), max_freq_ridges, max_chirp_ridges, ...
    fre_half, chirp_len, time_points, chirp_grid, win_len, fs, ...
    pcct_params, max_chirp_ridges_observed, 2);
end

function pcct_params = default_pcct_params()
pcct_params = struct();
pcct_params.init_points = 30;
pcct_params.freq_range = 5;
pcct_params.chirp_range = 5;
pcct_params.use_local_predict_score = true;
pcct_params.use_projection_amp_comp = true;
pcct_params.comp_global_energy_quantile = 0;
end

function params = fill_pcct_defaults(params)
if nargin < 1 || isempty(params)
    params = default_pcct_params();
end

defaults = default_pcct_params();
default_fields = fieldnames(defaults);

for i = 1:numel(default_fields)
    field_name = default_fields{i};
    if ~isfield(params, field_name) || isempty(params.(field_name))
        params.(field_name) = defaults.(field_name);
    end
end
end

function score = local_predict_score(patch, f_set, c_set, anchor_f, anchor_c, ...
    dt_local, freq_grid, chirp_grid, freq_range, chirp_range, global_tfcr_peak)
patch = double(patch);
patch_max = max(patch(:));
score = patch;

if patch_max <= 0
    return;
end

% Skip prediction steering when the local response is weak.
abs_gate = 0.08;
if patch_max < abs_gate * global_tfcr_peak
    return;
end

E = patch / (patch_max + eps);

f_anchor = freq_grid(anchor_f);
c_anchor = chirp_grid(anchor_c);
f_pred = f_anchor + 0.5 * c_anchor * dt_local;
c_pred = c_anchor;

F_grid = repmat(freq_grid(f_set).', 1, numel(c_set));
C_grid = repmat(chirp_grid(c_set), numel(f_set), 1);

df = abs(freq_grid(2) - freq_grid(1));
dc = abs(chirp_grid(2) - chirp_grid(1));

sigma_f = max(freq_range * df, eps);
sigma_c = max(chirp_range * dc, eps);

D_f = ((F_grid - f_pred) ./ sigma_f).^2;
D_c = ((C_grid - c_pred) ./ sigma_c).^2;

W = exp(-0.5 * D_f - 0.5 * D_c);

score = E;

strong_mask = E >= 0.75;
score(strong_mask) = E(strong_mask) .* ...
    (1.0 + 0.5 * W(strong_mask));
end

function TFR = render_tfr_from_ridges( ...
    tensor_mag, max_freq_ridges, max_chirp_ridges, ...
    fre_half, chirp_len, time_points, chirp_grid, win_len, fs, ...
    pcct_params, observed_chirp_ridges, spread_width)
if nargin < 11
    observed_chirp_ridges = max_chirp_ridges;
end
if nargin < 12 || isempty(spread_width)
    spread_width = 2;
end

num_modes = numel(max_freq_ridges);
TFR = zeros(time_points, fre_half);

use_amp_comp = pcct_params.use_projection_amp_comp;
sigma_t = (win_len / fs) / 6;
comp_global_energy_quantile = pcct_params.comp_global_energy_quantile;
amp_comp_energy_threshold = 0;

if use_amp_comp && comp_global_energy_quantile > 0
    ridge_energy_values = zeros(1, num_modes * time_points);
    ridge_energy_count = 0;

    for mode_idx = 1:num_modes
        f_trace = max_freq_ridges{mode_idx};
        c_trace = max_chirp_ridges{mode_idx};

        if isempty(f_trace) || isempty(c_trace)
            continue;
        end

        for nt = 1:time_points
            f_idx = max(1, min(fre_half, round(f_trace(nt))));
            c_idx = max(1, min(chirp_len, round(c_trace(nt))));
            ridge_energy_count = ridge_energy_count + 1;
            ridge_energy_values(ridge_energy_count) = tensor_mag(nt, f_idx, c_idx)^2;
        end
    end

    ridge_energy_values = ridge_energy_values(1:ridge_energy_count);
    ridge_energy_values = sort(ridge_energy_values(ridge_energy_values > 0));

    if ~isempty(ridge_energy_values)
        quantile_value = max(0, min(1, comp_global_energy_quantile));
        quantile_idx = max(1, ceil(quantile_value * numel(ridge_energy_values)));
        amp_comp_energy_threshold = ridge_energy_values(quantile_idx);
    end
end

for mode_idx = 1:num_modes
    f_trace = max_freq_ridges{mode_idx};
    c_trace = max_chirp_ridges{mode_idx};
    if mode_idx <= numel(observed_chirp_ridges)
        c_trace_observed = observed_chirp_ridges{mode_idx};
    else
        c_trace_observed = c_trace;
    end

    if isempty(f_trace) || isempty(c_trace)
        continue;
    end

    for nt = 1:time_points
        f_idx = max(1, min(fre_half, round(f_trace(nt))));
        c_idx = max(1, min(chirp_len, round(c_trace(nt))));
        c_idx_observed = max(1, min(chirp_len, round(c_trace_observed(nt))));
        amp_comp = 1;
        ridge_energy = tensor_mag(nt, f_idx, c_idx)^2;
        if use_amp_comp && ridge_energy >= amp_comp_energy_threshold
            delta_c = chirp_grid(c_idx_observed) - chirp_grid(c_idx);
            B_delta_c = 1 + 4 * pi^2 * sigma_t^4 * delta_c^2;
            amp_comp = B_delta_c^(1 / 4);
        end

        corrected_amp = tensor_mag(nt, f_idx, c_idx) * amp_comp;
        TFR = add_frequency_spread(TFR, nt, f_idx, corrected_amp, ...
            fre_half, spread_width);
    end
end
end

function ridge_s = smooth_ridge_indices(ridge, smooth_win, max_idx)
ridge = double(ridge);
if numel(ridge) > smooth_win
    if exist('smoothdata', 'file') == 2
        ridge_s = smoothdata(ridge, 'movmedian', smooth_win);
    else
        ridge_s = local_movmedian(ridge, smooth_win);
    end
else
    ridge_s = ridge;
end

ridge_s = round(ridge_s);
ridge_s = max(1, min(max_idx, ridge_s));
end

function y = local_movmedian(x, win_len)
y = x;
half_win = floor(win_len / 2);
for idx = 1:numel(x)
    left_idx = max(1, idx - half_win);
    right_idx = min(numel(x), idx + half_win);
    y(idx) = median(x(left_idx:right_idx));
end
end

function TFR = add_frequency_spread(TFR, nt, f_idx, amp_value, fre_half, spread_width)
for offset = -spread_width:spread_width
    current_f_idx = f_idx + offset;
    if current_f_idx < 1 || current_f_idx > fre_half
        continue;
    end

    if spread_width == 0
        weight = 1;
    else
        weight = exp(-(offset^2) / (2 * (spread_width / 2)^2));
    end

    TFR(nt, current_f_idx) = TFR(nt, current_f_idx) + weight * amp_value;
end
end

function z = analytic_signal_fft(x)
x = x(:).';
N = length(x);
X = fft(x);
h = zeros(1, N);

if mod(N, 2) == 0
    h(1) = 1;
    h(N / 2 + 1) = 1;
    h(2:N / 2) = 2;
else
    h(1) = 1;
    h(2:(N + 1) / 2) = 2;
end

z = ifft(X .* h);
end

