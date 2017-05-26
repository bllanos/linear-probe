function [ fg ] = plotHueClassifier( h_inc, h_classifiers, n_classes, varargin )
% PLOTHUECLASSIFIER  Visualize hue classifiers
%
% ## Syntax
% plotHueClassifier(...
%   h_inc, h_classifiers, n_classes [, legend_names, line_specs, fg]...
% )
% fg = plotHueClassifier(...
%   h_inc, h_classifiers, n_classes [, legend_names, line_specs, fg]...
% )
%
% ## Description
% plotHueClassifier(...
%   h_inc, h_classifiers, n_classes [, legend_names, line_specs, fg]...
% )
%   Creates a figure or updates an existing figure with plots of hue
%   classifiers, against a background of hue values.
%
% fg = plotHueClassifier(...
%   h_inc, h_classifiers, n_classes [, legend_names, line_specs, fg]...
% )
%   Additionally returns the handle of the figure which was created or updated.
%
% ## Input Arguments
%
% h_inc -- Increment between hue sample values
%   The hue spacing corresponding to adjacent values in the hue estimators.
%
% h_classifiers -- Evaluated hue density estimators
%   A horizontal concatenation of one or more instances of the `classifier`
%   output argument of `mlDiscreteClassifier`. In other words,
%   `h_classifiers(:,i)` contains the i-th classifier to plot.
%
% n_classes -- Number of classes
%   The number of hue classes, not including the "background class". The
%   background class is assumed to have an index of zero, and the classes
%   of interest have indices from one, up to, and including, `n_classes`.
%   `n_classes` determines the scale of the plot, and makes it obvious if
%   some classes are not used by the classifiers.
%
% legend_names -- Plot legend entries
%   A cell vector with a length equal to the number of columns in
%   `h_classifiers`, where the i-th cell contains the string to use as
%   the legend entry for the plot of `h_classifiers(:,i)`.
%
%   If empty or not passed, no legend is added to the plot.
%
% line_specs -- Plot line specifications
%   A cell vector where each cell contains a valid line spec string (i.e. a
%   string that could be the value of a 'LineSpec' argument passed to
%   `plot`).
%
%   Hue classifiers will be plotted with line specs that cycle through the
%   elements of `line_specs`.
%
%   Defaults to `{'k-'}` if empty or not passed.
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
%
% See also hueVariableKernelDensityEstimator, hueGaussianDensityEstimator, hsv2rgb, plot

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created May 26, 2017

nargoutchk(0, 1);
narginchk(3, 6);

if ~isempty(varargin)
    legend_names = varargin{1};
else
    legend_names = {};
end
if length(varargin) > 1
    line_specs = varargin{2};
else
    line_specs = {'k-'};
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

h = (0:h_inc:1).';

% Create a background rainbow
background = cat(3, repmat(h.', n_classes, 1), ones(n_classes, length(h), 2));
background = hsv2rgb(background);
ax = axes(fg);
image(ax, [0, 1], [0.5, n_classes], background)
ax.YDir = 'normal';

hold on
for i = 1:size(h_classifiers, 2)
    plot(...
            h, h_classifiers(:, i),...
            line_specs{mod(i - 1, length(line_specs)) + 1},...
            'LineWidth', 2.0...
        )
end
hold off

if ~isempty(legend_names)
    legend(legend_names{:});
end
xlabel('Hue, \theta (range [0, 1])')
ylabel('Class')

end

