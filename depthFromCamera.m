function [d] = depthFromCamera( P, X )
%DEPTHFROMCAMERA Depth of points in front of the principal plane
%
% ## Syntax
% [d] = depthFromCamera( P, X )
%
% ## Description
% [d] = depthFromCamera( P, X )
%   Returns the depth of points in front of the principal plane of a
%   finite projective camera.
%
% ## Input Arguments
%
% P -- Camera projection matrix
%   A 3 x 4 array representing a camera projection matrix. The camera must
%   be a finite projective camera, not an affine camera.
%
% X -- Point
%   An n x 3 array representing 3D points.
%
% ## Output Arguments
%
% d -- Depth
%   The distances of `X` from the center of the camera, measured
%   perpendicular to the camera's principal ray.

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created March 10, 2017

nargoutchk(1,1);
narginchk(2,2);

M = P(:, 1:3);
m3 = M(end, :);
X = [X ones(size(X, 1), 1)];
x = (P * X.').';
d = sign(det(M)) .* x(:, end) ./ (X(:, end) .* norm(m3));

end



