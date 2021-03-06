function [ H, S ] = rgb2hs( I )
% RGB2HS  Determine hue and saturation values from RGB values
%
% ## Syntax
% H = rgb2hs( I )
% [H, S] = rgb2hs( I )
%
% ## Description
% H = rgb2hs( I )
%   Returns the hue channel of the image.
% [H, S] = rgb2hs( I )
%   Additionally returns the saturation channel of the image.
%
% ## Input Arguments
%
% I -- Colour image
%   An m x n x 3 array representing an image in the RGB colour space.
%
% ## Output Arguments
%
% H -- Hue image
%   An m x n array representing the hue channel of `I`. Values are in the
%   range [0, 1].
%
% S -- Saturation image
%   An m x n array representing the saturation channel of `I`. Values are
%   in the range [0, 1].
%
% ## Notes
% - This is just a wrapper for the built-in `rgb2hsv` function. Originally,
%   I had written my own version (commented-out below), but verified that
%   its output differed by a very small amount (at most 1.12 hue degrees on
%   the image of a rainbow pattern) from the built-in function.
%
% ## References
% - T. Gevers and H. Stokman. "Robust Histogram Construction from Color
%   Invariants for Object Recognition". IEEE Transactions on Pattern
%   Analysis and Machine Intelligence, vol. 26, no. 1, pp. 113-118, Jan.
%   2004.
% - https://en.wikipedia.org/wiki/Hue
%
% See also rgb2hsv

% Bernard Llanos
% Spring 2016 research assistantship supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created August 9, 2016

nargoutchk(1, 2);
narginchk(1, 1);

% R = double(I(:, :, 1));
% G = double(I(:, :, 2));
% B = double(I(:, :, 3));
% H = atan2d(sqrt(3) * (G - B), 2 * R - G - B);
% H_neg_filter = (H < 0);
% H_neg = H(H_neg_filter);
% H(H_neg_filter) = 360 + H_neg;
% Output is in degrees, but can be normalized to the range [0, 1] by dividing by 360.

HS = rgb2hsv(I);
H = HS(:, :, 1);
if nargout > 1
    S = HS(:, :, 2);
end

end

