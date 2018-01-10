function [above, below] = reprojectProbe( lengths, widths, P, d, X_tip )
%REPROJECTPROBE Reproject measured points on the probe
%
% ## Syntax
% [above, below] = reprojectProbe( lengths, widths, P, d, X_tip )
%
% ## Description
% [above, below] = reprojectProbe( lengths, widths, P, d, X_tip )
%   Reproject the measured points on the probe onto the image plane.
%
% ## Input Arguments
%
% lengths -- Probe length measurements
%   A column vector of length `n` containing the measured physical
%   distances, from the tip of the probe, of points along the probe axis.
%
% widths -- Probe diameter measurements
%   A column vector of length `n` containing the measured diameters of the
%   probe at each of the locations corresponding to the points in
%   `lengths`.
%
% P -- Camera projection matrix
%   The camera matrix (intrinsic and extrinsic parameters). A 3 x 4 matrix.
%
% d -- Estimated probe axis
%   A unit row 3-vector containing the estimated direction of the probe in
%   space (pointing from the probe tip towards the other end of the probe).
%
% X_tip -- Probe tip
%   A row 3-vector containing the estimated position of the probe tip in
%   space.
%
% ## Output Arguments
%
% above -- Points along lower edge of probe
%   Points located at higher image y-coordinates than the probe axis in the
%   image. An n x 2 array of image coordinates of points on the contour of
%   the probe, at distances along the axis of the probe corresponding to
%   the elements of `lengths`.
%
% below -- Points along upper edge of probe
%   Points located at lower y-coordinates than the probe axis in the
%   image. An n x 2 array of image coordinates of points on the contour of
%   the probe, at distances along the axis of the probe corresponding to
%   the elements of `lengths`.
%
% See also planeNormalFromImageLine

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created March 27, 2017

nargoutchk(1, 2);
narginchk(5, 5);

d_image = (P * [d 0].').';
X_tip_image = (P * [X_tip 1].').';
image_line = cross(d_image, X_tip_image);
u = planeNormalFromImageLine(P, image_line);

n = size(lengths, 1);
X_tip = repmat(X_tip, n, 1);
u = repmat(u, n, 1);
axis = X_tip + repmat(lengths, 1, 3) .* repmat(d, n, 1);

r = repmat(widths / 2, 1, 3); % Take radii, not diameters
above = (P * [(axis + r .* u).'; ones(1, n)]).';
above = above(:, 1:2) ./ repmat(above(:, 3), 1, 2);
if nargout > 1
    below = (P * [(axis - r .* u).'; ones(1, n)]).';
    below = below(:, 1:2) ./ repmat(below(:, 3), 1, 2);
end

end

