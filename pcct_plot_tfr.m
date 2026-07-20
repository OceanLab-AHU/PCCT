function pcct_plot_tfr(tt, freq, tfr_data, color_limits, x_limits, y_limits, ...
    x_label, y_label, desired_gap, red_box, show_colorbar)
%PCCT_PLOT_TFR Shared PCCT time-frequency visualization.

if nargin < 9 || isempty(desired_gap)
    desired_gap = 0.172;
end
if nargin < 10
    red_box = [];
end
if nargin < 11 || isempty(show_colorbar)
    show_colorbar = true;
end

figure('Color', [1 1 1], 'Units', 'pixels', ...
    'Position', [20 100 338 230]);

[tfr_display, display_limits, use_unit_ticks] = ...
    prepare_tfr_display(tfr_data, color_limits);

[tt_plot, freq_plot, tfr_display] = interpolate_tfr_display( ...
    tt, freq, tfr_display, 2);
[tt_pc, freq_pc, tfr_pc] = pad_pcolor_grid(tt_plot, freq_plot, tfr_display);
[tt_mesh, freq_mesh] = meshgrid(tt_pc, freq_pc);
pcolor(tt_mesh, freq_mesh, tfr_pc);
shading interp;
caxis(display_limits);
colormap(pcct_white_colormap());

ax = gca;
xlabel(ax, x_label, 'FontName', 'Times New Roman', 'FontSize', 12);
ylabel(ax, y_label, 'FontName', 'Times New Roman', 'FontSize', 12);
set(ax, 'FontSize', 12, 'FontName', 'Times New Roman', ...
    'LineWidth', 1.5, 'GridLineStyle', 'none', 'Box', 'on', ...
    'TickLabelInterpreter', 'tex', 'YDir', 'normal');

if ~isempty(x_limits)
    xlim(x_limits);
end

if ~isempty(y_limits)
    ylim(y_limits);
end

if ~isempty(red_box)
    hold(ax, 'on');
    rectangle('Parent', ax, 'Position', red_box, ...
        'EdgeColor', [1 0 0], 'LineWidth', 1.4);
    hold(ax, 'off');
end

max_freq = max(abs(freq));
if max_freq >= 1e3
    ax.YAxis.Exponent = 3 * floor(log10(max_freq) / 3);
else
    ax.YAxis.Exponent = 0;
end

if show_colorbar
    c = colorbar(ax);
    c.Label.String = '';
    set(c, 'FontSize', 12, 'FontName', 'Times New Roman');
    if use_unit_ticks
        set(c, 'Ticks', 0:0.2:1);
    end

    ax_pos = get(ax, 'Position');
    c_width = 0.01;
    c.Position = [ax_pos(1) + ax_pos(3) + desired_gap, ...
        ax_pos(2), c_width, ax_pos(4)];
end
end

function [tt_plot, freq_plot, tfr_interp] = interpolate_tfr_display( ...
    tt, freq, tfr_display, interp_factor)
tt = tt(:).';
freq = freq(:).';

if numel(tt) < 2 || numel(freq) < 2 || interp_factor <= 1
    tt_plot = tt;
    freq_plot = freq;
    tfr_interp = tfr_display;
    return;
end

[X, Y] = meshgrid(tt, freq);
tt_plot = linspace(min(tt), max(tt), numel(tt) * interp_factor);
freq_plot = linspace(min(freq), max(freq), numel(freq) * interp_factor);
[Xq, Yq] = meshgrid(tt_plot, freq_plot);
tfr_interp = interp2(X, Y, tfr_display, Xq, Yq, 'linear', 0);
tfr_interp(~isfinite(tfr_interp)) = 0;
end

function [x_pad, y_pad, z_pad] = pad_pcolor_grid(x, y, z)
x = x(:).';
y = y(:).';
z_pad = z;

if numel(x) > 1
    dx = x(end) - x(end-1);
else
    dx = 1;
end
if numel(y) > 1
    dy = y(end) - y(end-1);
else
    dy = 1;
end

x_pad = [x, x(end) + dx];
y_pad = [y, y(end) + dy];
z_pad = [z_pad, z_pad(:, end)];
z_pad = [z_pad; z_pad(end, :)];
end

function [tfr_display, display_limits, use_unit_ticks] = ...
    prepare_tfr_display(tfr_data, color_limits)
tfr_display = tfr_data;
use_unit_ticks = false;

if isempty(color_limits)
    display_limits = [min(tfr_data(:)), max(tfr_data(:))];
    if display_limits(1) == display_limits(2)
        display_limits(2) = display_limits(1) + eps;
    end
    return;
end

display_limits = color_limits;
if numel(color_limits) == 2 && color_limits(1) >= 0 && color_limits(2) > 0
    tfr_display = tfr_data / color_limits(2);
    tfr_display(tfr_display < 0) = 0;
    tfr_display(tfr_display > 1) = 1;
    display_limits = [0 1];
    use_unit_ticks = true;
end
end

function cmap = pcct_white_colormap()
cmap = [
    ones(3, 3);
    linspace(1, 0.6, 10)', linspace(1, 0.8, 10)', ones(10, 1);
    linspace(0.6, 0, 10)', linspace(0.8, 0.3, 10)', ones(10, 1);
    zeros(10, 1), linspace(0.3, 0.8, 10)', ones(10, 1);
    zeros(10, 1), linspace(0.8, 1, 10)', linspace(1, 0.8, 10)';
    linspace(0, 1, 10)', ones(10, 1), linspace(0.8, 0, 10)';
    ones(10, 1), linspace(1, 0.5, 10)', zeros(10, 1);
    ones(10, 1), linspace(0.5, 0, 10)', zeros(10, 1);
    linspace(1, 0.5, 10)', zeros(10, 2)
    ];
end
