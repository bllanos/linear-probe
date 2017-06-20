function [ dist ] = hueGaussianDensityEstimator( H, mask, resolution )
% HUEGAUSSIANDENSITYESTIMATOR  Gaussian density estimator for hue values
%
% ## Syntax
% dist = hueGaussianDensityEstimator( H, mask, resolution )
%
% ## Description
% dist = hueGaussianDensityEstimator( H, mask, resolution )
%   Returns an evaluation of the Gaussian density estimator for hue values
%   at the given resolution.
%
% ## Input Arguments
%
% H -- Image hue channel
%   An h x w array, containing the hue values of the image (range [0, 1]).
%
% mask -- Region of interest
%   An h x w logical array, defining the region in the image for which the
%   Gaussian density estimator is to be computed.
%
% resolution -- Number of samples
%   The number of equally-spaced samples in the range of hue values from 0
%   (inclusive) to 1 (inclusive) at which to evaluate the Gaussian density
%   estimator.
%
%   Specifically, the evaluation is performed at `linspace(0,1,resolution)`.
%
% ## Output Arguments
%
% dist -- Evaluated Gaussian density estimator
%   A column vector of length equal to `resolution` containing the Gaussian
%   density estimator for the hue values in `H` evaluated at samples in the
%   interval [0, 1] as determined by `resolution`.
%
%   To find the approximate value of the estimator at a query value 'x' in
%   the range [0, 1], use `queryDiscretized1DFunction`.
%
% ## Notes
% - While a Gaussian density estimator is a continuous function,
%   computing a discrete version, by sampling, creates a lookup table for
%   fast approximate evaluation.
%
% ## References
% - T. Gevers and H. Stokman. "Robust Histogram Construction from Color
%   Invariants for Object Recognition". IEEE Transactions on Pattern
%   Analysis and Machine Intelligence, vol. 26, no. 1, pp. 113-118, Jan.
%   2004.
%
% See also rgb2hue, queryDiscretized1DFunction

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created May 26, 2017

nargoutchk(1, 1);
narginchk(3, 3);

h = H(mask);
mu = mean(h);
s = std(h);

x = linspace(0, 1, resolution);
% Theo Gevers and Harro Stokman use the following:
%   `deviation = mod(x(i) - h, 0.5);`
% This seems incorrect, and results in two peaks in the distribution.
deviation = abs(x - mu);
filter = deviation > 0.5;
deviation(filter) = 1.0 - deviation(filter);
dist = normpdf(deviation, 0, s);
normalization_factor = normcdf([-0.5 0.5], 0, s);
normalization_factor = diff(normalization_factor, 1, 2);
dist = dist / normalization_factor;

end

