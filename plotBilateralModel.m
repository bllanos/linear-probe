function [ fg ] = plotBilateralModel( model, varargin )
% PLOTBILATERALMODEL  Visualize the output of `bilateralModel`
%
% ## Syntax
% plotBilateralModel( model [, model_axes, image_size, colors, fg] )
% fg = plotBilateralModel( model [, model_axes, image_size, colors, fg] )
%
% ## Description
% plotBilateralModel( model [, model_axes, image_size, colors, fg] )
%   Creates a figure or updates an existing figure with the points in the
%   input structure created by `bilateralModel`.
%
% fg = plotBilateralModel( model [, model_axes, image_size, colors, fg] )
%   Additionally returns the handle of the figure which was created or updated.
%
% ## Input Arguments
%
% model -- Output of `bilateralModel`
%   The `model` or `model_px` output argument of `bilateralModel`, which is
%   a structure containing arrays of point coordinates.
%
% model_axes -- Lines representing PCA components
%   The `axes` output argument of `bilateralModel`, which contains the
%   coefficients of lines to be plotted.
%
%   If empty or not passed, no lines are plotted.
%
% image_size -- Size of image on which to draw lines representing PCA components
%   A two-element vector containing the pixel height and width of the image
%   on which the lines representing PCA components are to be drawn. The
%   figure is assumed to display an image, with the origin of the image at
%   the top-left corner of the figure.
%
%   When `model_axes` is passed and non-empty, `image_size` must be
%   non-empty.
%
%   Internally, `image_size` is passed to `lineToBorderPoints`.
%
% colors -- Plot colours
%   A cell vector of colours to use when plotting `model` and, if passed as
%   input, `model_axes`. Colours can be string colour codes recognized by
%   MATLAB, or numeric RGB colour value 3-vectors.
%
%   The first five elements specify plot colours for the points in the
%   fields of `model` in the following order:
%   - 'head'
%   - 'tail'
%   - 'above'
%   - 'below'
%   - 'unmatched'
%
%   The sixth and seventh elements, required only if `model_axes` is passed
%   and non-empty, are the colours to use for lines corresponding to the
%   first and second rows of `model_axes`, respectively.
%
%   If empty or not passed, `colors` defaults to
%   `{ 'b', 'k', 'g', 'r', 'y', 'c', 'm' }`.
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
% - The caller is responsible for adding axis labels, a plot title, and a
%   legend to the figure.
%
% See also bilateralModel, lineToBorderPoints

% Bernard Llanos
% Spring 2016 research assistantship supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created August 24, 2016

nargoutchk(0, 1);
narginchk(1, 5);

if ~isempty(varargin)
    model_axes = varargin{1};
else
    model_axes = [];
end
if length(varargin) > 1
    image_size = varargin{2};
else
    image_size = [];
end
if ~isempty(model_axes) && isempty(image_size)
    error('`image_size` is required when `model_axes` is non-empty.');
end
if length(varargin) > 2
    colors = varargin{3};
else
    colors = {};
end
if isempty(colors)
    colors = { 'b', 'k', 'g', 'r', 'y', 'c', 'm' };
end
if length(varargin) > 3
    fg = varargin{4};
    if ~ishandle(fg)
        error('The `fg` optional input argument must be a valid figure handle.');
    end
    figure(fg);
else
    fg = figure;
end

hold on
    
% Plot points
if isfield(model, 'head')
    head = model.head;
    scatter(head(1), head(2), colors{1})
end
if isfield(model, 'tail')
    tail = model.tail;
    scatter(tail(1), tail(2), colors{2})
end
above = model.above;
scatter(above(:, 1), above(:, 2), colors{3})
below = model.below;
scatter(below(:, 1), below(:, 2), colors{4})
unmatched = model.unmatched;
scatter(unmatched(:, 1), unmatched(:, 2), colors{5})

% Plot PCA lines
if ~isempty(model_axes)
    line_points = lineToBorderPoints(model_axes, image_size);
    line(line_points(1, [1,3])', line_points(1, [2,4])', 'Color', colors{6});
    line(line_points(2, [1,3])', line_points(2, [2,4])', 'Color', colors{7});
end

hold off

end

