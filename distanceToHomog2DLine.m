function [ d, s ] = distanceToHomog2DLine(p, l)
% DISTANCETOHOMOG2DLINE Perpendicular distance from points to a 2D homogenous line
%
% ## Syntax
% d = distanceToHomog2DLine(p, l)
% [ d, s ] = distanceToHomog2DLine(p, l)
%
% ## Description
% d = distanceToLine(p, l)
%   Returns the perpendicular distances from the points in `p` to the line.
%
% [ d, s ] = distanceToLine(p, l)
%   Additionally returns the sign of the distances.
%
% ## Input Arguments
%
% p -- Points
%   An n x 2 array of 2D points.
%
% l -- Line
%   A 1 x 3 array containing the homogenous representation of a 2D line.
%
% ## Output Arguments
%
% d -- Distances
%   An n x 1 array where `d(i)` is the perpendicular distance from `p(i,:)`
%   to the line.
%
% s -- Signs of distances
%   An n x 1 array, where `s(i)` represents which side of the line the
%   point `p(i, :)` is on. `s(i)` is the sign of the dot product of `p(i, :)`
%   and the normalized version of `l`.
%
% ## References
% - http://homepages.inf.ed.ac.uk/rbf/CVonline/LOCAL_COPIES/BEARDSLEY/node2.html
%   ("Manipulating Points and Lines," by Bob Fisher,
%    Fri Nov 7 12:08:26 GMT 1997)
%
% See also distanceToLine, closestPointOnLine

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created April 6, 2017

nargoutchk(1, 2);
narginchk(2, 2);

n = size(p, 1);
scale = sqrt(sum(l(1:2) .^ 2));
l = l ./ repmat(scale, 1, 3);
l = repmat(l, n, 1);

d = dot(l, [p ones(n, 1)], 2);

if nargout > 1
    s = sign(d);
end

d = abs(d);

end

