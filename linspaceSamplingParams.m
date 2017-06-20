function [ increment, resolution ] = linspaceSamplingParams( x, varargin )
% LINSPACESAMPLINGPARAMS  Retrieve the parameters used to sample a range
%
% ## Syntax
% increment = linspaceSamplingParams( x )
% increment = linspaceSamplingParams( a, b, resolution )
% [ increment, resolution ] = linspaceSamplingParams(____)
%
% ## Description
% increment = linspaceSamplingParams( x )
%   Returns the spacing between samples in `x`
%
% increment = linspaceSamplingParams( a, b, resolution )
%   Returns the spacing between samples implied by `a`, `b`, and
%   `resolution`
%
% [ increment, resolution ] = linspaceSamplingParams(____)
%   Additionally returns the number of samples
%
% ## Explanation
%
% This function recovers the missing information from the following
% function call: `x = linspace(a, b, resolution)`
%
% `x` is the output of a uniform sampling of a range, with a spacing of
% `increment` between consecutive samples (i.e. `all(diff(x) == increment)`
% is true).
%
% Given `x`, `resolution` and `increment` can be recovered. Given `a`, `b`,
% and `resolution`, `increment` can be recovered.
%
% See also linspace, queryDiscretized1DFunction

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created June 19, 2017

nargoutchk(1, 2);
narginchk(1, 3);

if isempty(varargin)
    a = x(1);
    b = x(end);
    resolution = length(x);
elseif length(varargin) == 2
    a = x;
    b = varargin{1};
    resolution = varargin{2};
else
    error('Unexpected number of input arguments.')
end

increment = (b - a) / (resolution - 1);

end