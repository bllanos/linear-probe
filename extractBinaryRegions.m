function [ regions ] = extractBinaryRegions( backprojected_distributions, radius, varargin )
% EXTRACTBINARYREGIONS  Threshold and obtain binary regions after filtering out poorly-distinguished pixels
%
% ## Syntax
% regions = extractBinaryRegions( backprojected_distributions, radius [, verbose] )
%
% ## Description
% regions = extractBinaryRegions( backprojected_distributions, radius [, verbose] )
%   Returns binary regions corresponding to each input backprojected
%   probability distribution.
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
% radius -- Radius of disk structuring element
%   The radius passed in `strel('disk', radius)` when creating a
%   structuring element to erode the final binary images prior to extraction
%   of connected components. Erosion is effective for removing noise, which
%   otherwise results in a large number of small connected components.
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
% ## Algorithm
% 
% To extract the regions which appear to correspond to a given colour
% class, as dictated by a probability distribution, it is sufficient to
% threshold the greyscale image formed by backprojecting the probability
% distribution. However, in this case, there are multiple colour classes,
% which are assumed to be quite distinct, such that pixels that have
% similar probabilities with respect to multiple colours actually do not
% belong to any colour. Consequently, such pixels are filtered out so that
% they are not taken into account when computing a histogram for Otsu
% thresholding, and are excluded from the binary connected components
% obtained after thresholding. These poorly-distinguished pixels are
% identified by Otsu thresholding applied to the pairwise absolute
% differences of the individual images in `backprojected_distributions`.
%
% See also bwconncomp, otsuthresh, imbinarize, imerode, ratioDistribution, hueVariableKernelDensityEstimator

% Bernard Llanos
% Spring 2016 research assistantship supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created August 15, 2016

nargoutchk(1, 1);
narginchk(2, 3);

if ~isempty(varargin)
    verbose = varargin{1};
else
    verbose = false;
end

image_height = size(backprojected_distributions, 1);
image_width = size(backprojected_distributions, 2);
n = size(backprojected_distributions, 3);

% For each colour, remove pixels which have the same values in the
% probability distributions computed for other colours.
pair_differentiation_masks = true(image_height, image_width, n);
for i = 1:n
    mask_i = pair_differentiation_masks(:, :, i);
    for j = (i + 1):n
        mask_j = pair_differentiation_masks(:, :, j);
        I_diff = imabsdiff(backprojected_distributions(:, :, i), backprojected_distributions(:, :, j));
        mask_ij = imbinarize(I_diff); % Otsu's method
        pair_differentiation_masks(:, :, i) = mask_i & mask_ij;
        pair_differentiation_masks(:, :, j) = mask_j & mask_ij;
    end    
end

% For the pixels that remain, use Otsu's method to obtain binary regions.
regions = struct('Connectivity', cell(n, 1), 'ImageSize', cell(n, 1),...
    'NumObjects', cell(n, 1), 'PixelIdxList', cell(n, 1));
for i = 1:n
    mask_i = pair_differentiation_masks(:, :, i);
    backprojected_distributions_i = backprojected_distributions(:, :, i);
    counts_i = imhist(backprojected_distributions_i(mask_i));
    threshold_i = otsuthresh(counts_i);
    bw_i = imbinarize(backprojected_distributions_i, threshold_i) & mask_i;
    disk = strel('disk', radius);
    bw_i = imerode(bw_i, disk);
    if verbose
        figure
        imshow(bw_i);
        title(sprintf('Binary image obtained for colour class %d', i))
    end
    regions(i) = bwconncomp(bw_i);
end

end

