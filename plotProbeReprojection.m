function [fg] = plotProbeReprojection( I, above, below, lengths, widths, P, d, X_tip, str )
%PLOTPROBEREPROJECTION Plot detected and reprojected points
%
% ## Syntax
% plotProbeReprojection(...
%     I, above, below, lengths, widths, P, d, X_tip, title...
% )
%
% fg = plotProbeReprojection(...
%     I, above, below, lengths, widths, P, d, X_tip, title...
% )
%
% ## Description
% plotProbeReprojection(...
%     I, above, below, lengths, widths, P, d, X_tip, title...
% )
%   Plot the points on the probe detected in the image and plot the
%   reprojections of the measured points on the probe.
%
% fg = plotProbeReprojection(...
%     I, above, below, lengths, widths, P, d, X_tip, title...
% )
%   Additionally return the figure.
%
% ## Input Arguments
%
% I -- Image
%   An image showing the probe, to allow for visual inspection of the
%   accuracy of probe localization. Points will be plotted on top of this
%   image.
%
% above -- Interest points along lower edge of probe
%   Points located at higher y-coordinates than the probe midline in the
%   image. An n x 2 array of image coordinates.
%
% below -- Interest points along upper edge of probe
%   Points located at lower y-coordinates than the probe midline in the
%   image. An n x 2 array of image coordinates. `below(i, :)` is the point
%   oppposite `above(i, :)` on an edge between two colour bands of the
%   probe.
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
%   The title to give to the figure.
%
% ## Output Arguments
%
% fg -- Figure handle
%   The handle to the figure containing the output.
%
% See also reprojectProbe, planeNormalFromImageLine

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created March 10, 2017

nargoutchk(0, 1);
narginchk(9, 9);

fg = figure;
imshow(I);

hold on

image_size = size(I);
image_size = image_size(1:2);
d_image = (P * [d 0].').';
X_tip_image = (P * [X_tip 1].').';
X_tip_image = X_tip_image ./ X_tip_image(end);
image_line = cross(d_image, X_tip_image);
line_points_plotting = lineToBorderPoints(image_line, image_size);
line(line_points_plotting([1,3]), line_points_plotting([2,4]), 'Color', 'c');

u = planeNormalFromImageLine(P, image_line);
u_image = (P * [u 0].').';
u_image_line = cross(u_image, X_tip_image);
line_points_plotting = lineToBorderPoints(u_image_line, image_size);
line(line_points_plotting([1,3]), line_points_plotting([2,4]), 'Color', 'r');

allPoints = [above; below];
scatter(allPoints(:, 1), allPoints(:, 2), 'g.');

[above, below, ~] = reprojectProbe( lengths, widths, P, d, X_tip );
allPoints = [above; below];
scatter(allPoints(:, 1), allPoints(:, 2), 'r.');

scatter(X_tip_image(1), X_tip_image(2), 'mo');

hold off
legend('Axis', 'Normal', 'Detected points', 'Reprojected 3D points', 'Tip');
title(str)

end

