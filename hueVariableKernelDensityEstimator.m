function [ dist, inc ] = hueVariableKernelDensityEstimator( H, R, G, B, mask, std_rgb, resolution )
% HUEVARIABLEKERNELDENSITYESTIMATOR  Variable kernel density estimator for hue values from RGB values
%
% ## Syntax
% dist = hueVariableKernelDensityEstimator(...
%   H, R, G, B, mask, std_rgb, resolution...
% )
% [ dist, inc ] = hueVariableKernelDensityEstimator(...
%   H, R, G, B, mask, std_rgb, resolution...
% )
%
% ## Description
% dist = hueVariableKernelDensityEstimator(...
%   H, R, G, B, mask, std_rgb, resolution...
% )
%   Returns an evaluation of the variable kernel density estimator for hue
%   values at the given resolution.
%
% [ dist, inc ] = hueVariableKernelDensityEstimator(...
%   H, R, G, B, mask, std_rgb, resolution...
% )
%   Returns an evaluation of the variable kernel density estimator for hue
%   values at the given resolution, and the spacing at which the estimator
%   was sampled.
%
% ## Input Arguments
%
% H -- Image hue channel
%   An h x w array, containing the hue values of the image (range [0, 1]).
%
% R -- Image red channel
%   An h x w array, containing the red values of the image (range [0, 1]).
%
% G -- Image green channel
%   An h x w array, containing the green values of the image (range [0, 1]).
%
% B -- Image blue channel
%   An h x w array, containing the blue values of the image (range [0, 1]).
%
% mask -- Region of interest
%   An h x w logical array, defining the region in the image for which the
%   variable kernel density estimator is to be computed.
%
% std_rgb -- RGB standard deviation fitted curves
%   Polynomial coefficients describing the relationship between RGB channel
%   values and RGB channel standard deviations.
%
%   This argument is of the form of the 'rgb_sigma_polyfit' variable output
%   by the script '.\EstimateRGBStandardDeviations.m'.
%
% resolution -- Number of samples
%   The number of equally-spaced samples in the range of hue values from 0
%   (inclusive) to 1 (inclusive) at which to evaluate the variable kernel
%   density estimator.
%
%   Specifically, the evaluation is performed at `linspace(0,1,resolution)`.
%
% ## Output Arguments
%
% dist -- Evaluated variable kernel density estimator
%   A column vector of length equal to `resolution` containing the variable
%   kernel density estimator for the hue values in `H` evaluated at samples
%   in the interval [0, 1] as determined by `resolution`.
%
%   The RGB channels of the image and their associated standard deviation
%   functions are used to modulate the scale parameter of the estimator.
%
% inc -- Increment between samples
%   The spacing between the samples in the range [0, 1] at which the
%   variable kernel density estimator has been evaluated.
%
%   To find the approximate value of the estimator at a query value 'x' in
%   the range [0, 1], use `queryDiscretized1DFunction`.
%
% ## Notes
% - While a variable kernel density estimator is a continuous function,
%   computing a discrete approximation to it by sampling has the following
%   advantages:
%   - The approximation can be computed in advance and then used as a
%     lookup table to reduce the number of operations in the future.
%   - The space required for storing the approximation can be much less
%     than that required to store all the input points used to compute the
%     exact (to within numerical error) evaluation of the estimator at a
%     query value.
% - It would be straight-forward to modify this function to accept an array
%   instead of the scalar `resolution` input argument, and evaluate the
%   variable kernel density estimator at the values in the array.
%
% ## References
% - T. Gevers and H. Stokman. "Robust Histogram Construction from Color
%   Invariants for Object Recognition". IEEE Transactions on Pattern
%   Analysis and Machine Intelligence, vol. 26, no. 1, pp. 113-118, Jan.
%   2004.
%
% See also rgb2hue, queryDiscretized1DFunction

% Bernard Llanos
% Spring 2016 research assistantship supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created August 9, 2016

    function [ s ] = hueSigma(r, g, b, std_rgb)
        std_r2 = polyval(std_rgb(:, 1), r) .^ 2;
        std_g2 = polyval(std_rgb(:, 2), g) .^ 2;
        std_b2 = polyval(std_rgb(:, 3), b) .^ 2;
        num = 3 * (...
                (b - g).^2 .* std_r2 +...
                (b - r).^2 .* std_g2 +...
                (r - g).^2 .* std_b2...
            );
        denom = 4 * (...
                r.^2 + g.^2 + b.^2 - r.*g - r.*b - g.*b...
            ) .^ 2;
        s = sqrt(num ./ denom);
        
        % Enforce a maximum value for numerical stability
        s(s > 0.5) = 0.5;
        s(isnan(s)) = 0.5; % Assume worst case for NaN values
        
        % Values of zero actually correspond to worst-case values as well;
        % The result of zero as opposed to NaN is due to numerical error
        s(s == 0) = 0.5;
    end

nargoutchk(2, 2);
narginchk(7, 7);

r = R(mask);
g = G(mask);
b = B(mask);
h = H(mask);
s = hueSigma(r, g, b, std_rgb);

x = linspace(0, 1, resolution);
dist = zeros(resolution, 1);
inc = 1 / (resolution - 1);
for i = 1:resolution
    % Theo Gevers and Harro Stokman use the following:
    %   `deviation = mod(x(i) - h, 0.5);`
    % This seems incorrect, and results in two peaks in the distribution
    % for each `h`.
    deviation = abs(x(i) - h);
    filter = deviation > 0.5;
    deviation(filter) = 1.0 - deviation(filter);
    p = normpdf(deviation, 0, s);
    dist(i) = sum(p);
end
dist = dist / length(h);

end

