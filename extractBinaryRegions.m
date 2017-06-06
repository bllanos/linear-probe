function [ regions, bw ] = extractBinaryRegions(...
    I, color_distributions, color_distribution_resolution,...
    rgb_sigma_polyfit, threshold, radius, varargin...
)
% EXTRACTBINARYREGIONS  Threshold and refine binary regions
%
% ## Syntax
% regions = extractBinaryRegions(...
%   I, color_distributions, color_distribution_resolution,...
%   rgb_sigma_polyfit, threshold, radius [, mask, verbose]...
% )
% [ regions, bw ] = extractBinaryRegions(...
%   I, color_distributions, color_distribution_resolution,...
%   rgb_sigma_polyfit, threshold, radius [, mask, verbose]...
% )
%
% ## Description
% regions = extractBinaryRegions(...
%   I, color_distributions, color_distribution_resolution,...
%   rgb_sigma_polyfit, threshold, radius [, mask, verbose]...
% )
%   Returns binary regions corresponding to each colour probability
%   distribution.
%
% [ regions, bw ] = extractBinaryRegions(...
%   I, color_distributions, color_distribution_resolution,...
%   rgb_sigma_polyfit, threshold, radius [, mask, verbose]...
% )
%   Additionally returns binary images corresponding to the binary regions.
%
% ## Input Arguments
%
% I -- Image
%   An RGB image in which the colours described by `color_distributions`
%   are to be detected.
%
% color_distributions -- Colour estimators
%   Discretized density estimators of image hue values corresponding to the
%   different colour classes. The i-th column of this 2D array stores the
%   estimator for the i-th colour class of to be detected in the image `I`.
%
% color_distribution_resolution -- Colour estimator sample count
%   The number of equally-spaced samples in the range of hue values from 0
%   (inclusive) to 1 (inclusive) at which the density estimators for colour
%   classes have been evaluated.
%
% rgb_sigma_polyfit -- Camera RGB noise model
%   An array describing the variation in RGB channel standard deviations
%   with RGB values in the image. This information should be computed from
%   images taken under the same conditions and with the same camera
%   parameters as the image (`I`), if not computed from this same image.
%
%   Refer to the documentation of './EstimateRGBStandardDeviations.m' for
%   details.
%
% threshold -- Threshold identifying poorly-distinguished pixels
%   The threshold applied to the probabilities assigned to pixels in `I`
%   from the distributions in `color_distributions` which eliminates pixels
%   which are not strongly selected by any distribution.
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
%   Defaults to `true(image_height,image_width)` if empty, or not passed.
%
% verbose -- Debugging flags
%   If recognized fields of `verbose` are true, corresponding graphical
%   output will be generated for debugging purposes.
%
%   All debugging flags default to false if `verbose` is not passed.
%
% ## Output Arguments
%
% regions -- Regions obtained by thresholding
%   A structure vector of length n, where the i-th element contains the
%   regions corresponding to the i-th colour class. Each element is of the
%   form of the output argument of `bwconncomp`.
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
% For this application, there are multiple colour classes. Therefore, the
% thresholding is performed on the colour class which gives the given pixel
% the highest probability. This is a simple approximation to a
% multi-labelling optimization problem (which could be solved with
% alpha-expansion, for example).
%
% See also bwconncomp, otsuthresh, imbinarize, imerode, mlDiscreteClassifier

% Bernard Llanos
% Spring 2016 research assistantship supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created August 15, 2016

nargoutchk(1, 2);
narginchk(6, 8);

image_width = size(I, 2);
image_height = size(I, 1);

mask = [];
if ~isempty(varargin)
    mask = varargin{1};
end
if isempty(mask)
    mask = true(image_height, image_width);
end
if length(varargin) > 1
    verbose = varargin{2};
    display_hue_image = verbose.display_hue_image;
    plot_hue_estimator = verbose.plot_hue_estimator;
    display_distribution_backprojections = verbose.display_distribution_backprojections;
    display_binary_images = verbose.display_binary_images;
else
    display_hue_image = false;
    plot_hue_estimator = false;
    display_distribution_backprojections = false;
    display_binary_images = false;
end

% Obtain hue values
H = rgb2hue(I);

I_double = im2double(I);
R = I_double(:, :, 1);
G = I_double(:, :, 2);
B = I_double(:, :, 3);

if display_hue_image
    figure
    H_color = ones(image_height, image_width, 3);
    H_color(:, :, 1) = H;
    H_color = hsv2rgb(H_color);
    imshowpair(H, H_color, 'montage');
    title('Hue channel of image')
end

% Compute the hue variable kernel density estimator for the image
[...
    I_color_distribution,...
    I_color_distribution_increment...
] = hueVariableKernelDensityEstimator(...
    H, R, G, B, mask,...
    rgb_sigma_polyfit, color_distribution_resolution...
);

n_colors = size(color_distributions, 2);

if plot_hue_estimator
    legend_names = cell(n_colors + 1, 1);
    legend_names{1} = 'Background';
    for i = 1:n_colors
        legend_names{i + 1} = sprintf('Probe colour %d', i);
    end
    plotHueDensityEstimator(...
        I_color_distribution_increment,...
        [I_color_distribution, color_distributions], legend_names...
    );
    title('Hue density estimators')
end

% Transform the image using histogram backprojection
n_colors_plus_background = n_colors + 1;
distributions_backprojected = zeros(image_height, image_width, n_colors_plus_background);
for i = 1:n_colors
    distributions_backprojected(:, :, i) = queryDiscretized1DFunction(...
            H, color_distributions(:, i), I_color_distribution_increment...
        );
end
distributions_backprojected(:, :, n_colors_plus_background) = queryDiscretized1DFunction(...
        H, I_color_distribution, I_color_distribution_increment...
    );

if display_distribution_backprojections
    for i = 1:n_colors
        figure
        imshow(distributions_backprojected(:, :, i));
        title(sprintf('Distribution backprojection for probe colour %d', i))
    end
    figure
    imshow(...
        distributions_backprojected(:, :, n_colors_plus_background) /...
        max(max(distributions_backprojected(:, :, n_colors_plus_background)))...
        );
    title('Distribution backprojection for the background')
end

use_otsu = isempty(threshold);

% Find the most likely colour class for each pixel
[~, max_ind] = max(distributions_backprojected, [], 3);
x = 0:(image_width - 1);
y = 1:image_height;
[X,Y] = meshgrid(x,y);
max_ind_linear = Y + X .* image_height + (max_ind - 1) .* (image_width * image_height);
bw = false(image_height, image_width, n_colors_plus_background);
bw(max_ind_linear) = true;
bw = bw(:, :, 1:(end - 1));

% Threshold probabilities for each colour
for i = 1:n_colors
    distributions_backprojected_i = distributions_backprojected(:, :, i);
    if use_otsu
        counts_i = imhist(distributions_backprojected_i(mask));
        threshold_i = otsuthresh(counts_i);
        bw(:, :, i) = bw(:, :, i) & imbinarize(distributions_backprojected_i, threshold_i);
    else
        bw(:, :, i) = bw(:, :, i) & imbinarize(distributions_backprojected_i, threshold);
    end
    % Apply the mask
    bw(:, :, i) = bw(:, :, i) & mask;
end

% Obtain binary regions.
regions = struct('Connectivity', cell(n_colors, 1), 'ImageSize', cell(n_colors, 1),...
    'NumObjects', cell(n_colors, 1), 'PixelIdxList', cell(n_colors, 1));
for i = 1:n_colors
    bw_i = bw(:, :, i);
    disk = strel('disk', radius);
    bw_i = imerode(bw_i, disk);
    if display_binary_images
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

