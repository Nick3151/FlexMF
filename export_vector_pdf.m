function export_vector_pdf(filename, fig)
%EXPORT_VECTOR_PDF  Save a figure as a true vector PDF.
%   export_vector_pdf(filename)
%   export_vector_pdf(filename, fig)
%
% Uses exportgraphics(..., 'ContentType','vector') so thin events from
% patch/line plots survive PDF export. Prefer this over save2pdf for
% plot_MR, SimpleWHPlot*_patch data panels, and other sparse rasters.
% save2pdf uses print -rDPI and can drop 1-bin features.

if nargin < 2 || isempty(fig)
    fig = gcf;
end
filename = char(filename);
if ~endsWith(filename, '.pdf', 'IgnoreCase', true)
    filename = [filename, '.pdf'];
end
if isgraphics(fig)
    set(fig, 'Renderer', 'painters');
end
exportgraphics(fig, filename, 'ContentType', 'vector');
end
