function [ mask ] = detectBoundingBoxes(...
    I, probe,...
    probe_color_distributions, rgb_sigma_polyfit,...
    params, varargin...
    )
%DETECTBOUNDINGBOXES Find a bounding box containing the probe
%
% ## Syntax
% mask = detectBoundingBoxes(...
%     I, probe,...
%     probe_color_distributions, rgb_sigma_polyfit,...
%     params [, verbose]...
% )
%
% ## Description
% mask = detectBoundingBoxes(...
%     I, probe,...
%     probe_color_distributions, rgb_sigma_polyfit,...
%     params [, verbose]...
% )
%   Returns a binary image describing the probe's bounding region.
%
% ## Input Arguments
%
% I -- Image containing the probe
%   An RGB image showing enough of the probe for detection to be possible.
%   The image has ideally been taken by the same camera, and under the same
%   camera parameters, as the image used to create the probe detection
%   model. Ideally, the lighting conditions in the image are similar to
%   those in the context of which the detection model was created.
%
%   The image should have been corrected for lens distortion.
%
% probe -- Probe measurements
%   Refer to the documentation of './CreateProbeDetectionModel.m' for
%   details.
%
% probe_color_distributions -- Probe colour estimators
%   Discretized density estimators of image hue values corresponding to the
%   different coloured bands on the probe, in the same order (starting from
%   the active tip of the probe). The i-th column of this 2D array stores
%   the estimator for the i-th colour class of probe segments.
%
% rgb_sigma_polyfit -- Camera RGB noise model
%   A parameter passed directly to 'extractBinaryRegions()'. As described
%   in the documentation of 'extractBinaryRegions()', an empty array can be
%   passed instead, to assume a uniform background colour distribution.
%
% params -- Fixed parameters
%   Parameters that should be stable across a variety of input images and
%   probe models. `params` is a structure containing the following fields:
%   - saturation_threshold: A minimum threshold on colour saturation
%       values, passed as the `saturation_threshold` parameter of
%       'extractBinaryRegions()'
%   - erosion_radius: Radius for morphological erosion of the binary images
%       describing the probe colour regions. The `radius` parameter of
%       'extractBinaryRegions()'
%   - radius_adj: Pixel radius used to filter candidate probe colour
%       regions to those close to regions for other colours. The
%       `radius_adj` parameter of 'detectProbeBinaryRegions()'.
%   - axis_distance_outlier_threshold: Number of standard deviations from
%       the estimate of the probe axis beyond which a region is determined
%       to be distinct from the probe. The
%       `axis_distance_outlier_threshold` parameter of
%       'detectProbeBinaryRegions()'.
%   - dilation_radius: Radius for dilating the candidate probe colour
%       regions to include some of the background. This parameter is used
%       as part of bounding region construction.
%   - region_expansion_factor_length: Factor by which to expand oriented
%       boxes fitted to the dilated candidate probe colour regions, along
%       their major axis directions.
%   - region_expansion_factor_width: Factor by which to expand oriented
%       boxes fitted to the dilated candidate probe colour regions, along
%       their minor axis directions.
%
% verbose -- Debugging flags
%   If recognized fields of `verbose` are true, corresponding graphical
%   output will be generated for debugging purposes.
%
%   All debugging flags default to false if `verbose` is not passed.
%
% ## Output Arguments
%
% mask -- Probe bounding region
%   A binary image (with the same width and height as `I`), where nonzero
%   pixels indicate bounding areas for the coloured bands of the probe. In
%   general, there are multiple, disconnected regions.
%
% See also hueVariableKernelDensityEstimator, hueGaussianDensityEstimator, extractBinaryRegions, detectProbeBinaryRegions

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created March 28, 2017

%% Parse input arguments

nargoutchk(1,1);
narginchk(5,6);

if ~isempty(varargin)
    verbose = varargin{1};
    display_original_image = verbose.display_original_image;
    extractBinaryRegionsVerbose.display_hue_image = verbose.display_hue_image;
    extractBinaryRegionsVerbose.display_saturation_image = verbose.display_saturation_image;
    extractBinaryRegionsVerbose.plot_hue_estimator = verbose.plot_hue_estimator;
    extractBinaryRegionsVerbose.plot_hue_classifier = verbose.plot_hue_classifier;
    extractBinaryRegionsVerbose.display_distribution_backprojections = verbose.display_distribution_backprojections;
    extractBinaryRegionsVerbose.display_binary_images = verbose.display_binary_images;
    verbose_region_filtering = verbose.verbose_region_filtering;
    display_region_expansion = verbose.display_region_expansion;
else
    display_original_image = false;
    extractBinaryRegionsVerbose.display_hue_image = false;
    extractBinaryRegionsVerbose.display_saturation_image = false;
    extractBinaryRegionsVerbose.plot_hue_estimator = false;
    extractBinaryRegionsVerbose.plot_hue_classifier = false;
    extractBinaryRegionsVerbose.display_distribution_backprojections = false;
    extractBinaryRegionsVerbose.display_binary_images = false;
    verbose_region_filtering = false;
    display_region_expansion = false;
end

image_width = size(I, 2);
image_height = size(I, 1);
image_n_channels = size(I, 3);
if image_n_channels ~= 3
    error('A 3-channel RGB image is required.')
end
if display_original_image
    figure;
    imshow(I);
    title('Base image');
end

%% Find an initial bounding area for the probe

[ probe_regions_initial, probe_regions_bw_initial] = extractBinaryRegions(...
    I, probe_color_distributions,...
    params.saturation_threshold,...
    rgb_sigma_polyfit,...
    params.erosion_radius,...
    [],...
    extractBinaryRegionsVerbose...
);
n_colors = size(probe_color_distributions, 2);
probe_regions_initial = probe_regions_initial(1:n_colors);
probe_regions_bw_initial = probe_regions_bw_initial(:, :, 1:n_colors);

% Identify which pairs of colours can be adjacent
probe_color_pairs = [probe.colors(1:(end - 1)), probe.colors(2:end)];
probe_color_pairs = probe_color_pairs(all(probe_color_pairs ~= 0, 2), :);
probe_color_pairs = sort(probe_color_pairs, 2);
probe_color_pairs = unique(probe_color_pairs, 'rows');

[
    ~,...
    probe_regions_bw_initial_filtered...
] = detectProbeBinaryRegions(...
        probe_regions_initial,...
        probe_regions_bw_initial,...
        probe_color_pairs,...
        params.radius_adj,...
        params.axis_distance_outlier_threshold,...
        verbose_region_filtering...
    );

% Obtain an upper bound on the probe's area
probe_regions_bw_initial_all = any(probe_regions_bw_initial_filtered, 3);
if display_region_expansion
    figure
    imshow(probe_regions_bw_initial_all);
    title('Initial probe colour regions prior to expansion')
end
disk = strel('disk', params.dilation_radius);
probe_regions_bw_initial_all = imdilate(probe_regions_bw_initial_all, disk);
% Find oriented bounding boxes
initial_region_bounds = regionprops(probe_regions_bw_initial_all, {...
    'Centroid', 'Orientation', 'MajorAxisLength', 'MinorAxisLength'...
    });
n_initial_regions = length(initial_region_bounds);
initial_region_rect_corners = cell(n_initial_regions, 1);
for i = 1:n_initial_regions
    half_length = initial_region_bounds(i).MajorAxisLength * params.region_expansion_factor_length / 2;
    half_width = initial_region_bounds(i).MinorAxisLength * params.region_expansion_factor_width / 2;
    angle = initial_region_bounds(i).Orientation;
    length_vector = half_length * [ cosd(angle), -sind(angle) ];
    width_vector = half_width * [ cosd(angle + 90), -sind(angle + 90) ];
    center = initial_region_bounds(i).Centroid;
    initial_region_rect_corners{i} = [
        center + length_vector + width_vector;
        center - length_vector + width_vector;
        center - length_vector - width_vector;
        center + length_vector - width_vector;
        ];
end

% Obtain a mask for the initial bounding boxes
mask = false(image_height, image_width);
for i = 1:n_initial_regions
    mask = mask | roipoly(...
            I, initial_region_rect_corners{i}(:, 1), initial_region_rect_corners{i}(:, 2)...
        );
end
if display_region_expansion
    figure
    probe_regions_bw_initial_all_display = repmat(probe_regions_bw_initial_all, 1, 1, 3);
    probe_regions_bw_initial_all_display(:, :, 1) =...
        probe_regions_bw_initial_all_display(:, :, 1) |...
        imdilate(bwperim(mask), strel('disk', 2));
    imshow(double(probe_regions_bw_initial_all_display));
    title('Initial probe colour regions after dilation, with outline of initial mask')
end
end

