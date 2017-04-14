function [d, u] = probeOrientationFromEndpoints( P, X_tip, X_end )
%PROBEORIENTATIONFROMENDPOINTS Compute probe orientation vectors given its endpoints
%
% ## Syntax
% [d, u] = probeOrientationFromEndpoints( length, P, X_tip, X_end )
%
% ## Description
% [d, u] = probeOrientationFromEndpoints( length, P, X_tip, X_end )
%   Calculate the axis and normal vectors of the probe, given the positions
%   of its endpoints.
%
% ## Input Arguments
%
% P -- Camera projection matrix
%   The camera matrix (intrinsic and extrinsic parameters).
%
% X_tip -- Probe tip
%   The 3D position of one point along the axis of the probe, as a 3 x 1
%   vector.
%
% X_tip -- Probe end
%   The 3D position of another point on axis of the probe, as a 3 x 1
%   vector.
%
% ## Output Arguments
%
% d -- Probe orientation vector
%   A 3 x 1 unit vector aligned with the axis of the probe in space.
%
% u -- Probe normal vector
%   A 3 x 1 unit vector perpendicular to the plane formed by the axis of
%   the probe in space and the camera centre, pointing towards the bottom
%   of the image.
%
% ## References
% - R. Hartley and A. Zisserman. Multiple View Geometry in Computer Vision,
%   2nd Edition. Cambridge, UK: Cambridge University Press, 2003.
%
% See also planeNormalFromImageLine

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created March 27, 2017

nargoutchk(1, 2);
narginchk(3, 3);

X_tip_image = (P * [X_tip 1].').';
X_end_image = (P * [X_end 1].').';
d = X_end - X_tip;
d = d ./ repmat(norm(d), 1, 3); % Normalize

image_line = cross(X_end_image, X_tip_image);
if nargout > 1
    u = planeNormalFromImageLine(P, image_line);
end

end