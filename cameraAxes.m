function [x, y, z] = cameraAxes( P )
%CAMERAAXES Find the axes of the image in world coordinates
%
% ## Syntax
% [x, y, z] = cameraAxes( P )
%
% ## Description
% [x, y, z] = cameraAxes( P )
%   Returns the directions in world coordinates corresponding to the
%   vanishing points of the image axes, and the camera's principal ray.
%
% ## Input Arguments
%
% P -- Camera projection matrix
%   A 3 x 4 array representing a camera projection matrix.
%
% ## Output Arguments
%
% x -- Image x-axis
%   The direction in world 3D coordinates corresponding to the image
%   positive x-axis. A 3 x 1 vector.
%
% y -- Image y-axis
%   The direction in world 3D coordinates corresponding to the image
%   positive y-axis. A 3 x 1 vector.
%
% z -- Principal ray
%   The direction in world 3D coordinates corresponding to the camera's
%   principal ray. A 3 x 1 vector.

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created March 10, 2017

nargoutchk(1,3);
narginchk(1,1);

M = P(:, 1:3);
m3 = M(end, :);
z = det(M) * m3;
z = z.' ./ norm(z);
P_inv = pinv(P); % Pseudoinverse
x = P_inv * [1; 0; 0];
x = x(1:3) ./ norm(x);
y = P_inv * [0; 1; 0];
y = y(1:3) ./ norm(y);

end

