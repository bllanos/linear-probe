function [u] = planeNormalFromImageLine(P, image_line)
%planeNormalFromImageLine Find the normal to a world plane given its line in the image
%
% ## Syntax
% u = planeNormalFromImageLine(image_line, P)
%
% ## Description
% u = planeNormalFromImageLine(image_line, P)
%   Determines the normal vector to the plane projected as the given line
%   in the image from the given camera.
%
% ## Input Arguments
%
% P -- Camera projection matrix
%   The camera matrix (intrinsic and extrinsic parameters)
%
% image_line -- Projection of plane
%   A 3-vector containing the homogenous coordinates of a line in the image
%   corresponding to the world plane.
%
% ## Output Arguments
%
% u -- Plane normal vector
%   A 3 x 1 unit vector perpendicular to the plane formed by the image line
%   and the camera centre, pointing towards the bottom of the image.
%
% ## References
% - R. Hartley and A. Zisserman. Multiple View Geometry in Computer Vision,
%   2nd Edition. Cambridge, UK: Cambridge University Press, 2003.

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created March 27, 2017

nargoutchk(1, 1);
narginchk(2, 2);

% The plane back-projected through the image line: Result 8.2 from Multiple
% View Geometry in Computer Vision, 2nd Edition
plane = image_line * P;

% Normal vector to the plane
u = plane(1:3);
u = u ./ repmat(norm(u), 1, 3); % Normalize
% Make sure the vector points to the bottom of the image
u_image = (P * [u 0].').';
u_image = u_image(1:2) ./ repmat(u_image(end), 1, 2);
y_centerline = -(image_line(1) * u_image(1) + image_line(3)) / image_line(2);
if y_centerline > u_image(2)
    u = -u;
end
end