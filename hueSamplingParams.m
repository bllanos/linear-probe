function [ increment, resolution ] = hueSamplingParams( x )
% HUESAMPLINGPARAMS  Retrieve the parameters used to sample the range of hues
%
% ## Syntax
% increment = hueSamplingParams( x )
% increment = hueSamplingParams( resolution )
% [ increment, resolution ] = hueSamplingParams(____)
%
% ## Description
% increment = hueSamplingParams( x )
%   Returns the spacing between hue samples used to generate `x`
%
% increment = hueSamplingParams( resolution )
%   Returns the spacing between hue samples implied by `resolution`
%
% [ increment, resolution ] = hueSamplingParams(____)
%   Additionally returns the number of hue samples.
%
% ## Explanation
%
% Given a vector `x` (of length greater than one), this function is
% equivalent to `[...] = linspaceSamplingParams(0, 1, length(x))`.
%
% Given a scalar `resolution`, this function is equivalent to
% `[...] = linspaceSamplingParams(0, 1, resolution)`.
%
% `x` is either a set of uniformly-spaced hue samples, or a function
% evaluated at a set of uniformly-spaced hue samples. The range of hue
% values is assumed to be from zero to one (inclusive).
%
% `resolution` is the number of uniformly-spaced hue samples.
%
% See also linspaceSamplingParams, queryDiscretized1DFunction

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created June 19, 2017

nargoutchk(1, 2);
narginchk(1, 1);

if length(x) > 1
    [increment, resolution] = linspaceSamplingParams(0, 1, length(x));
else
    [increment, resolution] = linspaceSamplingParams(0, 1, x);
end

end