function [fg] = plotProbeReprojection( I, above, below, lengths, widths, P, d, u, X_tip, str )
%PLOTPROBEREPROJECTION Linear estimation of probe tip and orientation
%
% ## Syntax
% plotProbeReprojection(...
%     I, above, below, lengths, widths, P, d, u, X_tip, title...
% )
%
% fg = plotProbeReprojection(...
%     I, above, below, lengths, widths, P, d, u, X_tip, title...
% )
%
% ## Description
% plotProbeReprojection(...
%     I, above, below, lengths, widths, P, d, u, X_tip, title...
% )
%   Plot the points on the probe detected in the image and plot the
%   reprojections of the measured points on the probe.
%
% fg = plotProbeReprojection(...
%     I, above, below, lengths, widths, P, d, u, X_tip, title...
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
%   the image points in `above` and `below`. The camera must be a finite
%   camera, not an affine camera.
%
% d -- Estimated probe axis
%   A unit 3-vector containing the estimated direction of the probe in
%   space (pointing from the probe tip towards the other end of the probe).
%
% u -- Estimated probe normal vector
%   A unit 3-vector containing the estimated direction in space relating
%   the pairs of points in `above` and `below`.
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

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created March 10, 2017

nargoutchk(0, 1);
narginchk(10, 10);

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

u_image = (P * [u 0].').';
u_image_line = cross(u_image, X_tip_image);
line_points_plotting = lineToBorderPoints(u_image_line, image_size);
line(line_points_plotting([1,3]), line_points_plotting([2,4]), 'Color', 'r');

n = size(above, 1);
nAll = 2 * n;
allPoints = [above; below];
scatter(allPoints(:, 1), allPoints(:, 2), 'g.');

l = repmat(lengths, 2, 1);
r = repmat(widths / 2, 2, 1); % Take radii, not diameters
r(n+1:end) = -r(n+1:end); % Account for the opposition between `above` and `below`
reprojected_points = (P * (repmat([X_tip 1], nAll, 1) + l .* repmat([d 0], nAll, 1) + r .* repmat([u 0], nAll, 1)).').';
reprojected_points = reprojected_points(:, 1:2) ./ repmat(reprojected_points(:, 3), 1, 2);
scatter(reprojected_points(:, 1), reprojected_points(:, 2), 'r.');

scatter(X_tip_image(1), X_tip_image(2), 'mo');

hold off
legend('Axis', 'Normal', 'Detected points', 'Reprojected 3D points', 'Tip');
title(str)

end

