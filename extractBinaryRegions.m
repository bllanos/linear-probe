function [ regions, bw ] = extractBinaryRegions( backprojected_distributions, threshold, radius, varargin )
% EXTRACTBINARYREGIONS  Threshold and refine binary regions
%
% ## Syntax
% regions = extractBinaryRegions(...
%   backprojected_distributions, threshold, radius [, mask, verbose]...
% )
% [ regions, bw ] = extractBinaryRegions(...
%   backprojected_distributions, threshold, radius [, mask, verbose]...
% )
%
% ## Description
% regions = extractBinaryRegions(...
%   backprojected_distributions, threshold, radius [, mask, verbose]...
% )
%   Returns binary regions corresponding to each input backprojected
%   probability distribution.
%
% [ regions, bw ] = extractBinaryRegions(...
%   backprojected_distributions, threshold, radius [, mask, verbose]...
% )
%   Additionally returns binary images corresponding to the binary regions.
%
% ## Input Arguments
%
% backprojected_distributions -- Backprojected probability distributions
%   An image_height x image_width x n array, where
%   `backprojected_distributions(:,:,i)` is a greyscale image where the
%   value of a pixel is drawn from the probability distribution of the i-th
%   colour. In other words, `backprojected_distributions(:,:,i)` is an
%   image where each pixel's value is the likelihood that the pixel belongs
%   to a region of the i-th colour.
%
% threshold -- Threshold identifying poorly-distinguished pixels
%   The threshold applied to the probabilities in
%   `backprojected_distributions` to eliminate pixels which are not
%   strongly selected by any distribution.
%
%   If empty, Otsu's method will be used to select the threshold
%   automatically.
%
%   Refer the discussion of the algorithm below for details.
%
% radius -- Radius of disk structuring element
%   The radius passed in `strel('disk', radius)` when creating a
%   structuring element to erode the final binary images prior to extraction
%   of connected components. Erosion is effective for removing noise, which
%   otherwise results in a large number of small connected components.
%
% mask -- Region of interest
%   A logical array of size image_height x image_width which determines the
%   set of pixels operated on by the function.
%
%   Defaults to `true(image_height,image_width)` if empty or not passed.
%
% verbose -- Debugging flag
%   If true, graphical output will be generated for debugging purposes.
%
%   Defaults to false if not passed.
%
% ## Output Arguments
%
% regions -- Regions obtained by thresholding
%   A structure vector of length n, where the i-th element contains the
%   regions corresponding to the i-th colour. Each element is of the form
%   of the output argument of `bwconncomp`.
%
% bw - Binary images
%   An image_height x image_width x n logical array, where `bw(:,:,i)` is
%   the binary image consisting of the connected components in
%   `regions(i)`.
%
% ## Algorithm
% 
% To extract the regions which appear to correspond to a given colour
% class, as dictated by a probability distribution, it is sufficient to
% threshold the greyscale image formed by backprojecting the probability
% distribution. This function does so, optionally choosing the threshold
% using Otsu's method.
%
% However, for this application, there are multiple colour classes.
% Therefore, the thresholding is performed on the colour class which gives
% the given pixel the highest probability. This is a simple approximation
% to a multi-labelling optimization problem (which could be solved with
% alpha-expansion, for example).
%
% See also bwconncomp, otsuthresh, imbinarize, imerode, mlDiscreteClassifier

% Bernard Llanos
% Spring 2016 research assistantship supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created August 15, 2016

nargoutchk(1, 2);
narginchk(3, 5);

image_height = size(backprojected_distributions, 1);
image_width = size(backprojected_distributions, 2);
n = size(backprojected_distributions, 3);

mask = [];
if ~isempty(varargin)
    mask = varargin{1};
end
if isempty(mask)
    mask = true(image_height, image_width);
end
if length(varargin) > 1
    verbose = varargin{2};
else
    verbose = false;
end

use_otsu = isempty(threshold);

% Find the most likely colour class for each pixel
[~, max_ind] = max(backprojected_distributions, [], 3);
x = 0:(image_width - 1);
y = 1:image_height;
[X,Y] = meshgrid(x,y);
max_ind_linear = Y + X .* image_height + (max_ind - 1) .* (image_width * image_height);
bw = false(image_height, image_width, n);
bw(max_ind_linear) = true;

% Threshold probabilities for each colour
for i = 1:n
    backprojected_distributions_i = backprojected_distributions(:, :, i);
    if use_otsu
        counts_i = imhist(backprojected_distributions_i(mask));
        threshold_i = otsuthresh(counts_i);
        bw(:, :, i) = bw(:, :, i) & imbinarize(backprojected_distributions_i, threshold_i);
    else
        bw(:, :, i) = bw(:, :, i) & imbinarize(backprojected_distributions_i, threshold);
    end
    % Apply the mask
    bw(:, :, i) = bw(:, :, i) & mask;
end

% Obtain binary regions.
regions = struct('Connectivity', cell(n, 1), 'ImageSize', cell(n, 1),...
    'NumObjects', cell(n, 1), 'PixelIdxList', cell(n, 1));
for i = 1:n
    bw_i = bw(:, :, i);
    disk = strel('disk', radius);
    bw_i = imerode(bw_i, disk);
    if verbose
        figure
        imshow(bw_i);
        title(sprintf('Binary image obtained for colour class %d', i))
    end
    regions(i) = bwconncomp(bw_i);
    if nargout > 1
        bw(:, :, i) = bw_i;
    end
end

end

