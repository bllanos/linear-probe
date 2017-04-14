function [ p_filter ] = isOnBothSidesOfLine( p, l )
% ISONBOTHSIDESOFLINE Test if regions cross a line
%
% ## Syntax
% p_filter = isOnBothSidesOfLine(p, l)
% p_filter = isOnBothSidesOfLine(p, line_points)
%
% ## Description
% p_filter = isOnBothSidesOfLine(p, l)
%   Returns a filter for regions on both sides of the line represented as a
%   homogenous vector.
%
% p_filter = isOnBothSidesOfLine(p, line_points)
%   Returns a filter for regions on both sides of the line represented by
%   two endpoints.
%
% ## Input Arguments
%
% p -- Test points
%   An cell vector of two dimensional arrays, where `p{i}` is an n_i x 2
%   array containing the coordinates of 2D points.
%
% l -- Line
%   A 1 x 3 array containing the homogenous representation of a 2D line.
%
% line_points -- Points defining the line
%   A 2 x 2 array containing the two points defining the line as its rows.
%
% ## Output Arguments
%
% p_filter -- Region filter
%   A logical vector, with the same length as `p`, where the i-th element
%   is true if there are points in `p{i}` on both sides of the line `l`.
%
% See also distanceToHomog2DLine, distanceToLine

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created April 6, 2017

    function b = isOnBothSides(signs)
        b = sum(signs < 0) & sum(signs >= 0);
    end

nargoutchk(1, 1);
narginchk(2, 2);

n = size(p, 1);

p_filter = false(n, 1);
if all(size(l) == [1 3])
    % Line defined by a homogenous equation
    for i = 1:n
        [~, signs] = distanceToHomog2DLine(p{i}, l);
        p_filter(i) = isOnBothSides(signs);
    end
elseif all(size(l) == [2 2])
    % Line defined by two endpoints
    for i = 1:n
        [~, signs] = distanceToLine(p{i}, l);
        p_filter(i) = isOnBothSides(signs);
    end
else
    error('Unknown form of the second input argument `l` or `line_points`.')
end
end

