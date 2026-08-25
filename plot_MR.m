function plot_MR(M, R, title_suffix)
%PLOT_MR  Display residual R (top) and motion field M (bottom).
%   plot_MR(M, R)
%   plot_MR(M, R, title_suffix)  appends title_suffix to subplot titles
%                                as 'R — …' and 'M — …'.

if nargin < 3 || isempty(title_suffix)
    title_R = 'R';
    title_M = 'M';
else
    title_R = ['R — ', char(title_suffix)];
    title_M = ['M — ', char(title_suffix)];
end

cmap_red = [ones(128,1),linspace(1,0,128)',linspace(1,0,128)'];
cmap_blue = [linspace(0,1,128)',linspace(0,1,128)',ones(128,1)];
cmap = [cmap_blue; cmap_red];
epsilon = 1e-4;

colormap(cmap)
ax_res = subplot('Position', [0.05, 0.55, 0.8, 0.4]);
maxValue = max(abs([R;M]), [], 'all')+epsilon;
clims = [-maxValue, maxValue]; 
imagesc(R, clims)
title(title_R, 'FontSize', 16, 'Interpreter', 'none')
set(ax_res, 'XTickLabel', [], 'YTickLabel', []);
colorbar('Position', [0.9 0.55 0.05 0.4], 'FontSize', 14);

ax_flux = subplot('Position', [0.05, 0.05, 0.8, 0.4]);
clims = [-maxValue, maxValue]; 
imagesc(M, clims)
title(title_M, 'FontSize', 16, 'Interpreter', 'none')
set(ax_flux, 'XTickLabel', [], 'YTickLabel', []);
colorbar('Position', [0.9 0.05 0.05 0.4], 'FontSize', 14);
set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])
