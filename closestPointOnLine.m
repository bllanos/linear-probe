function q = closestPointOnLine( l, p )
%CLOSESTPOINTONLINE Closest point on a line to the given point
%
% ## Syntax
% q = closestPointOnLine( l, p )
%
% ## Description
% q = closestPointOnLine( l, p )
%   Returns the closest points on lines `l` to the points `p`.
%
% ## Input Arguments
%
% l -- Lines
%   A n x 3 array containing the homogenous representations of 2D lines.
%
% p -- Points
%   The inhomogenous representation of 2D points; An n x 2 array.
%
% ## Output Arguments
%
% q -- Points on the lines
%   The closest points on the lines `l` to the points `p`. An n x 2 array.

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created March 9, 2017

nargoutchk(1,1);
narginchk(2,2);

% Derivation:
% p = sym('p', [2, 1]);
% l = sym('l', [3, 1]);
% n = l(1:2); % Normal vector to the line
% q = cross(cross([p; 1], [n; 0]), l)

q = [
    l(:,1).*l(:,3) + l(:,2).*(l(:,1).*p(:,2) - l(:,2).*p(:,1)),...
    l(:,2).*l(:,3) - l(:,1).*(l(:,1).*p(:,2) - l(:,2).*p(:,1)),...
    - l(:,1).^2 - l(:,2).^2
    ];
q = q(:, 1:2) ./ repmat(q(:, 3), 1, 2);

end

