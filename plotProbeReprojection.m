function [I_out] = plotProbeReprojection( I, points, lengths, widths, P, d, X_tip, varargin )
%PLOTPROBEREPROJECTION Plot detected and reprojected points
%
% ## Syntax
% plotProbeReprojection(...
%     I, points, lengths, widths, P, d, X_tip, title...
% )
%
% I_out = plotProbeReprojection(...
%     I, points, lengths, widths, P, d, X_tip...
% )
%
% ## Description
% plotProbeReprojection(...
%     I, points, lengths, widths, P, d, X_tip, title...
% )
%   Plot the points on the probe detected in the image and plot the
%   reprojections of the measured points on the probe.
%
% I_out = plotProbeReprojection(...
%     I, points, lengths, widths, P, d, X_tip...
% )
%   Produce an annotated image instead of plotting the information in a
%   figure. No legend or title are added to the image, in contrast to the
%   other invocation syntax.
%
% ## Input Arguments
%
% I -- Image
%   An image showing the probe, to allow for visual inspection of the
%   accuracy of probe localization. Points will be plotted on top of this
%   image.
%
% points -- Interest points along the edges of the probe
%   Points detected or annotated along the probe's contour in the
%   image. An n x 2 array of image coordinates.
%
%   An empty array can be passed.
%
% lengths -- Probe length measurements
%   A vector of length `n` containing the measured physical distances of
%   the pairs of points in `above` and `below` from the tip of the probe.
%
% widths -- Probe diameter measurements
%   A vector of length `n` containing the measured diameters of the probe
%   at each of the locations corresponding to the points in `above` and
%   `below`.
%
% P -- Camera projection matrix
%   The camera matrix (intrinsic and extrinsic parameters) corresponding to
%   the image points in `above` and `below`.
%
% d -- Estimated probe axis
%   A unit 3-vector containing the estimated direction of the probe in
%   space (pointing from the probe tip towards the other end of the probe).
%
% X_tip -- Probe tip
%   A 3-vector containing the estimated position of the probe tip in space.
%
% str -- Figure title
%   The title to give to the figure, if the output is to be a figure
%   instead of an annotated image.
%
% ## Output Arguments
%
% I_out -- Annotated image
%   A version of `I` with annotations added to show the location of the
%   probe.
%
%   The annotations on the image are as follows:
%   - cyan line: Estimated probe axis (`d`)
%   - red line: Line normal to the probe axis, passing through the probe
%       tip
%   - green circles: The locations in `points`
%   - red x's: Reprojected points on the edges between the bands of the
%       probe
%   - magenta star: Reprojection of the probe's tip
%
% ## Notes
% - When in image annotation mode (`I_out` is produced, instead of a
%   figure), points outside the image boundaries will not appear in the
%   output.
%
% See also reprojectProbe, planeNormalFromImageLine

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created March 10, 2017

nargoutchk(0, 1);
use_figure = (nargout == 0);
if use_figure
    narginchk(8, 8);
else
    narginchk(7, 7);
end

% Calculations required to produce annotations, or plot elements
image_size = size(I);
image_size = image_size(1:2);
d_image = (P * [d 0].').';
X_tip_image = (P * [X_tip 1].').';
X_tip_image = X_tip_image ./ X_tip_image(end);
image_line = cross(d_image, X_tip_image);
line_points_plotting_d = lineToBorderPoints(image_line, image_size);

u = planeNormalFromImageLine(P, image_line);
u_image = (P * [u 0].').';
u_image_line = cross(u_image, X_tip_image);
line_points_plotting_u = lineToBorderPoints(u_image_line, image_size);

have_detected_points = ~isempty(points);

[above, below] = reprojectProbe( lengths, widths, P, d, X_tip );
allPoints = [above; below];

% Output
if use_figure
    figure;
    imshow(I);

    hold on
    line(line_points_plotting_d([1,3]), line_points_plotting_d([2,4]), 'Color', 'c');
    line(line_points_plotting_u([1,3]), line_points_plotting_u([2,4]), 'Color', 'r');
    if have_detected_points
        scatter(points(:, 1), points(:, 2), 'g.');
    end
    scatter(allPoints(:, 1), allPoints(:, 2), 'r.');
    scatter(X_tip_image(1), X_tip_image(2), 'mo');
    hold off
    
    if have_detected_points
        legend('Axis', 'Normal', 'Detected points', 'Reprojected 3D points', 'Tip');
    else
        legend('Axis', 'Normal', 'Reprojected 3D points', 'Tip');
    end
    title(varargin{1})
else
    
    line_width = ceil(max(2, 0.002 * max(image_size)));
    marker_size = ceil(max(5, 0.01 * max(image_size)));
    
    I_out = insertShape(...
        I, 'line', line_points_plotting_d,...
        'LineWidth', line_width, 'Color', 'cyan'...
    );
    I_out = insertShape(...
        I_out, 'line', line_points_plotting_u,...
        'LineWidth', line_width, 'Color', 'red'...
    );
    if have_detected_points
        points_clipped = points(...
            points(:, 1) > 0 & points(:, 2) > 0 &...
            points(:, 1) <= image_size(2) & points(:, 2) < image_size(1),...
            :...
        );
        I_out = insertMarker(...
            I_out, points_clipped, 'o', 'Color', 'green', 'size', marker_size...
        );
    end
    allPoints_clipped = allPoints(...
        allPoints(:, 1) > 0 & allPoints(:, 2) > 0 &...
        allPoints(:, 1) <= image_size(2) & allPoints(:, 2) < image_size(1),...
        :...
    );
    I_out = insertMarker(...
        I_out, allPoints_clipped, 'x', 'Color', 'red', 'size', marker_size...
    );
    if X_tip_image(1) > 0 && X_tip_image(2) > 0 &&...
       X_tip_image(1) <= image_size(2) && X_tip_image(2) < image_size(1)
        I_out = insertMarker(...
            I_out, X_tip_image(1:2), '*', 'Color', 'magenta', 'size', marker_size...
        );
    end
end

end

