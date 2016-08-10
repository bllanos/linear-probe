function [ r ] = ratioDistribution( fg, bg )
% RATIODISTRIBUTION  Ratio distribution of signal and combined signal and background distributions
%
% ## Syntax
% r = ratioDistribution( fg, bg )
%
% ## Description
% r = ratioDistribution( fg, bg )
%   Returns the distribution describing the probability that a value is
%   part of the foreground signal given the distributions of the foreground
%   signal and of the entire dataset containing the foreground signal.
%
% ## Input Arguments
%
% fg -- Foreground distribution
%   A vector containing the probability distribution of values in the
%   foreground signal.
%
% bg -- Background and foreground distribution
%   A vector the same dimensions as `fg` containing the probability
%   distribution of values in the combined foreground and environment data.
%
% ## Output Arguments
%
% r -- Ratio distribution
%   A vector the same dimensions as `fg` containing the probabilities of
%   membership in the foreground signal for values in the combined
%   foreground and background given the foreground and combined foreground
%   and background distributions.
%
% ## Notes
% - The input distributions should be comparable, in that they should be
%   at the same scales and result from evaluating probability densities at
%   the same values of the sample variables.
%
% ## References
% - M.J. Swain and D.H. Ballard, "Color Indexing". International Journal of
%   Computer Vision, vol. 7, no. 1, pp. 11-32, 1991.

% Bernard Llanos
% Spring 2016 research assistantship supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created August 10, 2016

nargoutchk(1, 1);
narginchk(2, 2);

% Note that `min` will select `1` when `1` is paired with `Inf` and `NaN` values.
r = min(fg ./ bg, ones(size(fg)));
end

