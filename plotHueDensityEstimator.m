function [ fg ] = plotHueDensityEstimator( h_distributions, varargin )
% PLOTHUEDENSITYESTIMATOR  Visualize hue density estimators
%
% ## Syntax
% plotHueDensityEstimator(...
%   h_distributions [, legend_names, line_styles, fg]...
% )
% fg = plotHueDensityEstimator(...
%   h_distributions [, legend_names, line_styles, fg]...
% )
%
% ## Description
% plotHueDensityEstimator(...
%   h_distributions [, legend_names, line_styles, fg]...
% )
%   Creates a figure or updates an existing figure with plots of hue
%   density estimators.
%
% fg = plotHueDensityEstimator(...
%   h_distributions [, legend_names, line_styles, fg]...
% )
%   Additionally returns the handle of the figure which was created or updated.
%
% ## Input Arguments
%
% h_distributions -- Evaluated hue density estimators
%   A horizontal concatenation of one or more instances of the `dist`
%   output argument of `hueVariableKernelDensityEstimator` or
%   `hueGaussianDensityEstimator`, for example. In other words,
%   `h_distributions(:,i)` contains the i-th estimator to plot.
%
% legend_names -- Plot legend entries
%   A cell vector with a length equal to the number of columns in
%   `h_distributions`, where the i-th cell contains the string to use as
%   the legend entry for the plot of `h_distributions(:,i)`.
%
%   If empty or not passed, no legend is added to the plot.
%
% line_styles -- Plot line styles
%   A cell vector where each cell contains a valid line style string (i.e.
%   a string that could be the value of the 'LineStyle' name-value pair
%   argument passed to `plot`).
%
%   Hue estimators will be plotted with line styles that cycle through the
%   elements of `line_styles`.
%
%   Defaults to `{'-'}` if empty or not passed.
%
% fg -- Figure handle
%   A handle to the figure to update.
%
%   If not passed, a new figure is created.
%
% ## Output Arguments
%
% fg -- Figure handle
%   A handle to the figure updated with the plot output.
%
% ## Notes
% - The caller is responsible for adding a plot title. Axis labels are
%   added by this function.
% - The colour with which an estimator is plotted corresponds to the
%   maximum value of the estimator.
%
% See also hueVariableKernelDensityEstimator, hueGaussianDensityEstimator, hsv2rgb, plot

% Bernard Llanos
% Spring 2016 research assistantship supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created August 24, 2016

nargoutchk(0, 1);
narginchk(1, 4);

if ~isempty(varargin)
    legend_names = varargin{1};
else
    legend_names = {};
end
if length(varargin) > 1
    line_styles = varargin{2};
else
    line_styles = {'-'};
end
if length(varargin) > 2
    fg = varargin{3};
    if ~ishandle(fg)
        error('The `fg` optional input argument must be a valid figure handle.');
    end
    figure(fg);
else
    fg = figure;
end

h_inc = hueSamplingParams( h_distributions(:, 1) );
h = (0:h_inc:1).';

hold on
for i = 1:size(h_distributions, 2)
    h_distribution_i = h_distributions(:, i);
    h_distribution_i_normalized = h_distribution_i / sum(h_distribution_i);
    [~, color_ind] = max(h_distribution_i_normalized);
    color = h(color_ind);
    color = hsv2rgb([color, 1, 1]);
    plot(...
            h, h_distribution_i,...
            'Color', color,...
            'LineStyle', line_styles{mod(i - 1, length(line_styles)) + 1},...
            'LineWidth', 2.0...
        )
end
hold off

if ~isempty(legend_names)
    legend(legend_names{:});
end
xlabel('Hue, \theta (range [0, 1])')
ylabel('Density, P(\theta)')

end

