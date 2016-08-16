function [ d ] = distanceToLine(p, line_points)
% DISTANCETOLINE Perpendicular distance from points to a line defined by two points
%
% ## Syntax
% d = distanceToLine(p, line_points)
%
% ## Description
% d = distanceToLine(p, line_points)
%   Returns the perpendicular distances from the points in `p` to the line.
%
% ## Input Arguments
%
% p -- Test points
%   An n x k array of query points, where 'n' is the number of points, and
%   'k' is the dimensionality of each point.
%
% line_points -- Points defining the line
%   A 2 x k array containing the two points defining the line.
%
% ## Output Arguments
%
% d -- Distances
%   An n x 1 array where `d(i)` is the perpendicular distance from `p(i,:)`
%   to the line through `line_points`.
%
% ## Notes
% - The results of this function are unchanged if the rows of `line_points` are
%   exchanged.
%
% ## References
% - https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line

% Bernard Llanos
% Spring 2016 research assistantship supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created August 15, 2016

nargoutchk(1, 1);
narginchk(2, 2);

start = line_points(1, :);
endpoint = line_points(2, :);
direction = endpoint - start;
dimensionality = size(direction, 2);
unit_direction = direction ./ repmat(sqrt(dot(direction, direction, 2)), 1, dimensionality);

n = size(p, 1);
start_to_p = p - repmat(start, n, 1);
unit_direction = repmat(unit_direction, n, 1);
d = start_to_p - (repmat(dot(start_to_p, unit_direction, 2), 1, dimensionality) .* unit_direction);
d = sqrt(dot(d, d, 2));

end

