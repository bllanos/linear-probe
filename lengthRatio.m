function [ ratio ] = lengthRatio( points )
% LENGTHRATIO  Length ratio of three 1D points
%
% ## Syntax
% ratio = lengthRatio( points )
%
% ## Description
% ratio = lengthRatio( points )
%   Returns the length ratio of the three points.
%
% ## Input Arguments
%
% points -- Points on the affine line
%   A 3-vector, where `points(i)` contains the coordinate of the i-th 1D point.
%
% ## Output Arguments
%
% ratio -- Length ratio
%   A scalar value equal to the length ratio of `points`.
%
% ## Notes
% - The length ratio computed by this formula will change if the order
%   of the input points is reversed. (Contrast with the cross ratio of four
%   points.)
%
% See also crossRatio

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created March 29, 2017

ratio = (points(2) - points(1)) / (points(3) - points(2));

end

