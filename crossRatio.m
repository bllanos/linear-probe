function [ ratio ] = crossRatio( points )
% CROSSRATIO  Cross ratio of four 1D points
%
% ## Syntax
% ratio = crossRatio( points )
%
% ## Description
% ratio = crossRatio( points )
%   Returns the cross ratio of the four points.
%
% ## Input Arguments
%
% points -- Points on the projective line
%   A 4-vector, where `points(i)` contains the coordinate of the i-th 1D point.
%
% ## Output Arguments
%
% ratio -- Cross ratio
%   A scalar value equal to the cross ratio of `points`.
%
% ## Notes
% - The cross ratio computed by this formula does not change if the order
%   of the input points is reversed.
%
% ## References
% Page 45 of Hartley, Richard, and Andrew Zisserman. Multiple View Geometry In Computer Vision.
%   Cambridge, UK: Cambridge University Press, 2003. eBook Academic Collection
%   (EBSCOhost). Web. 27 July 2016.

% Bernard Llanos
% Spring 2016 research assistantship supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created July 27, 2016

ratio = ((points(1) - points(2)) * (points(3) - points(4)))...
    / ((points(1) - points(3)) * (points(2) - points(4)));

end

