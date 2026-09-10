function plot_MR(M, R, title_suffix, style)
%PLOT_MR  Display residual R (top) and motion field M (bottom).
%   plot_MR(M, R)
%   plot_MR(M, R, title_suffix)
%   plot_MR(M, R, title_suffix, style)
%
% title_suffix : appended to subplot titles as 'R — …' / 'M — …'
% style        : 'patch' (default) | 'imagesc'
%   'patch'   — per-neuron vector patches (good for FlexMF M/R on data X)
%   'imagesc' — classic heatmap with colorbar (good for EMD between two
%               short simulated sequences, e.g. EMD_demo / EMD_demo2)

if nargin < 3 || isempty(title_suffix)
    title_R = 'R';
    title_M = 'M';
else
    title_R = ['R — ', char(title_suffix)];
    title_M = ['M — ', char(title_suffix)];
end
if nargin < 4 || isempty(style)
    style = 'patch';
end
style = lower(char(style));
assert(ismember(style, {'patch', 'imagesc'}), ...
    'plot_MR: style must be ''patch'' or ''imagesc''.');

cmap_red = [ones(128,1),linspace(1,0,128)',linspace(1,0,128)'];
cmap_blue = [linspace(0,1,128)',linspace(0,1,128)',ones(128,1)];
cmap = [cmap_blue; cmap_red];
epsilon = 1e-4;
maxValue = max(abs([R; M]), [], 'all') + epsilon;
[N, T] = size(R);
clims = [-maxValue, maxValue];

switch style
    case 'imagesc'
        ax_res = subplot('Position', [0.05, 0.55, 0.8, 0.4]);
        imagesc(R, clims)
        title(title_R, 'FontSize', 16, 'Interpreter', 'none')
        set(ax_res, 'XTickLabel', [], 'YTickLabel', []);
        colormap(ax_res, cmap);
        colorbar('Position', [0.9 0.55 0.05 0.4], 'FontSize', 14);

        ax_flux = subplot('Position', [0.05, 0.05, 0.8, 0.4]);
        imagesc(M, clims)
        title(title_M, 'FontSize', 16, 'Interpreter', 'none')
        set(ax_flux, 'XTickLabel', [], 'YTickLabel', []);
        colormap(ax_flux, cmap);
        colorbar('Position', [0.9 0.05 0.05 0.4], 'FontSize', 14);

        set(gcf, 'Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8])

    case 'patch'
        ax_res = subplot('Position', [0.05, 0.55, 0.9, 0.4]);
        plot_signed_rows(ax_res, R, maxValue);
        title(title_R, 'FontSize', 16, 'Interpreter', 'none')
        set(ax_res, 'XTickLabel', [], 'YTickLabel', [], 'Box', 'on', ...
            'XLim', [0.5, T+0.5], 'YLim', [0, N+1], 'YDir', 'normal');

        ax_flux = subplot('Position', [0.05, 0.05, 0.9, 0.4]);
        plot_signed_rows(ax_flux, M, maxValue);
        title(title_M, 'FontSize', 16, 'Interpreter', 'none')
        set(ax_flux, 'XTickLabel', [], 'YTickLabel', [], 'Box', 'on', ...
            'XLim', [0.5, T+0.5], 'YLim', [0, N+1], 'YDir', 'normal');

        set(gcf, 'Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8], ...
            'Color', 'w', 'InvertHardcopy', 'off', 'Renderer', 'painters')
end
end

function plot_signed_rows(ax, Matrix, maxValue)
% Per-neuron vector patches: positive=red, negative=blue, both drawn upward
% with height |value| so pos/neg sit on the same vertical scale.
% Neuron 1 is drawn at the top (imagesc convention).
[N, T] = size(Matrix);
Xs = [1, 1:T, T];
dn = 1;
scaled = Matrix ./ maxValue;
cla(ax);
hold(ax, 'on');
for ni = 1:N
    base = dn * (N - ni + 1) - dn / 2;
    row = scaled(ni, :);
    row_pos = max(row, 0);
    row_neg = min(row, 0);
    patch(ax, Xs, [base, base + row_pos, base], [1 0 0], ...
        'EdgeColor', 'none', 'FaceAlpha', 1);
    % row_neg < 0: base - row_neg = base + |row_neg| (same upward scale as pos)
    patch(ax, Xs, [base, base - row_neg, base], [0 0 1], ...
        'EdgeColor', 'none', 'FaceAlpha', 1);
end
hold(ax, 'off');
end
