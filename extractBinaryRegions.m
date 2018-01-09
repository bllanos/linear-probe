function [ regions, bw ] = extractBinaryRegions(...
    I, color_distributions,...
    rgb_sigma_polyfit, radius, varargin...
)
% EXTRACTBINARYREGIONS  Threshold and refine binary regions
%
% ## Syntax
% regions = extractBinaryRegions(...
%   I, color_distributions,...
%   rgb_sigma_polyfit, radius [, mask, verbose]...
% )
% regions = extractBinaryRegions(...
%   I, color_distributions,...
%   [], radius [, mask, verbose]...
% )
% [ regions, bw ] = extractBinaryRegions(____)
%
% ## Description
% regions = extractBinaryRegions(...
%   I, color_distributions,...
%   rgb_sigma_polyfit, radius [, mask, verbose]...
% )
%   Returns binary regions corresponding to each colour probability
%   distribution, taking into account the background colour probability
%   distribution
%
% regions = extractBinaryRegions(...
%   I, color_distributions,...
%   [], radius [, mask, verbose]...
% )
%   Returns binary regions corresponding to each colour probability
%   distribution, assuming a uniform background colour probability
%   distribution
%
% [ regions, bw ] = extractBinaryRegions(____)
%   Additionally returns binary images corresponding to the binary regions.
%
% ## Input Arguments
%
% I -- Image
%   An RGB image, in which the colours described by `color_distributions`
%   are to be detected.
%
% color_distributions -- Colour estimators
%   Discretized density estimators of image hue values corresponding to the
%   different colour classes. The i-th column of this 2D array stores the
%   estimator for the i-th colour class of to be detected in the image `I`.
%
% rgb_sigma_polyfit -- Camera RGB noise model
%   An array describing the variation in RGB channel standard deviations
%   with RGB values in the image. This information should be computed from
%   images taken under the same conditions and with the same camera
%   parameters as the image (`I`), if not computed from this same image.
%
%   Refer to the documentation of 'noise_estimation/EstimateRGBStandardDeviations*.m' for
%   details.
%
%   `rgb_sigma_polyfit` is used to compute a probability distribution for
%   the colours in the image `I`. The resulting distribution is compared
%   with the colour estimators in `color_distributions` to improve colour
%   classification. If an empty array (`[]`) is passed instead of
%   `rgb_sigma_polyfit`, a uniform distribution will be used instead during
%   colour classification. Consequently, colour classification may be less
%   accurate, but computation time will be considerably lower.
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
% See also bwconncomp, imerode, mlDiscreteClassifier

% Bernard Llanos
% Spring 2016 research assistantship supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created August 15, 2016

nargoutchk(1, 2);
narginchk(5, 6);

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
    plot_hue_classifier = verbose.plot_hue_classifier;
    display_distribution_backprojections = verbose.display_distribution_backprojections;
    display_binary_images = verbose.display_binary_images;
else
    display_hue_image = false;
    plot_hue_estimator = false;
    plot_hue_classifier = false;
    display_distribution_backprojections = false;
    display_binary_images = false;
end

% Obtain hue values
H = rgb2hue(I);

if display_hue_image
    figure
    H_color = ones(image_height, image_width, 3);
    H_color(:, :, 1) = H;
    H_color = hsv2rgb(H_color);
    imshowpair(H, H_color, 'montage');
    title('Hue channel of image')
end

uniform_background = isempty(rgb_sigma_polyfit);

if ~uniform_background
    I_double = im2double(I);
    R = I_double(:, :, 1);
    G = I_double(:, :, 2);
    B = I_double(:, :, 3);

    % Compute the hue variable kernel density estimator for the image
    [...
        color_distribution_increment,...
        color_distribution_resolution...
    ] = hueSamplingParams( color_distributions(:, 1) );
    I_color_distribution = hueVariableKernelDensityEstimator(...
        H, R, G, B, mask,...
        rgb_sigma_polyfit, color_distribution_resolution...
    );
else
    [...
        color_distribution_increment,...
    ] = hueSamplingParams( color_distributions(:, 1) );
end

n_colors = size(color_distributions, 2);

if plot_hue_estimator
    legend_names = cell(n_colors + 1, 1);
    legend_names{1} = 'Background';
    for i = 1:n_colors
        legend_names{i + 1} = sprintf('Probe colour %d', i);
    end
    if uniform_background
        I_color_distribution = ones(size(color_distributions, 1), 1);
    end
    plotHueDensityEstimator(...
        [I_color_distribution, color_distributions], legend_names...
    );
    title('Hue density estimators')
end

if display_distribution_backprojections
    if uniform_background
        [color_classifier, color_classification_likelihood] = mlDiscreteClassifier(...
            color_distributions, color_distribution_increment, 'periodic'...
        );
    else
        [color_classifier, color_classification_likelihood] = mlDiscreteClassifier(...
            color_distributions, color_distribution_increment, I_color_distribution...
        );
    end
else
    if uniform_background
        color_classifier = mlDiscreteClassifier(...
            color_distributions, color_distribution_increment, 'periodic'...
        );
    else
        color_classifier = mlDiscreteClassifier(...
            color_distributions, color_distribution_increment, I_color_distribution...
        );
    end
end

if plot_hue_classifier
    plotHueClassifier(color_classifier, n_colors);
    title('Hue classifier for colors on the probe vs. the interest area')
end

% Find the most likely colour class for each pixel
H_masked = H(mask);
max_ind_masked = queryDiscretized1DFunction(...
    H_masked, color_classifier, color_distribution_increment...
);
max_ind = zeros(image_height, image_width);
max_ind(mask) = max_ind_masked;
foreground_mask = (max_ind ~= 0);
x = 0:(image_width - 1);
y = 1:image_height;
[X,Y] = meshgrid(x,y);
max_ind_linear = Y + X .* image_height + (max_ind - 1) .* (image_width * image_height);
max_ind_linear_masked = max_ind_linear(foreground_mask);
bw = false(image_height, image_width, n_colors);
bw(max_ind_linear_masked) = true;

% Transform the image using histogram backprojection
if display_distribution_backprojections
    max_probability_masked = queryDiscretized1DFunction(...
        H_masked, color_classification_likelihood, color_distribution_increment...
    );
    max_probability_masked = max_probability_masked(max_ind_masked ~= 0);
    distributions_backprojected = zeros(image_height, image_width, n_colors);
    distributions_backprojected(max_ind_linear_masked) = max_probability_masked;
    for i = 1:n_colors
        figure
        imshow(distributions_backprojected(:, :, i));
        title(sprintf('Distribution backprojection for probe colour %d', i))
    end
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

