%% Probe Object Detection
% Use previously estimated colour distributions and measured dimensions of
% the probe object to locate it in a test image.
%
% ## Usage
%   Modify the parameters and the paths to input data in the first code
%   section below, then run.
%
% ## Probe Detection Assumptions and Limitations
% - The image contains fairly saturated colours throughout, to assist with
%   hue-based colour discrimination. At minimum, unsaturated areas are far
%   from the probe.
% - The image only contains a single instance of the probe.
% - The probe may be occluded in one or more places. However, it can be
%   detected if at least 5 edges between pairs of coloured bands of the
%   probe are completely visible. The edges do not need to be between bands
%   that are all consecutive.
%
% ## Input
%
% ### Probe detection model
% A '.mat' file containing several variables, which is the output of
% '.\CreateProbeDetectionModel.m'. Refer to the documentation of
% '.\CreateProbeDetectionModel.m' for details.
%
% ### Image containing the probe
% An RGB image, in any format that can be loaded by `imread`, showing enough
% of the probe for detection to be possible. The image has ideally been
% taken by the same camera, and under the same camera parameters, as the
% image used to create the probe detection model. Ideally, the lighting
% conditions in the image are similar to those in the context of which the
% detection model was created.
%
% The image should have been corrected for lens distortion.
%
% ### Colour noise parameters
% A '.mat' file containing a 'rgb_sigma_polyfit' variable, as output by the
% script '.\EstimateRGBStandardDeviations.m'. 'rgb_sigma_polyfit' describes
% the variation in RGB channel standard deviations with RGB values in the
% image. This information should be computed from images taken
% under the same conditions and with the same camera parameters as the
% image in which the probe is to be detected, if not computed from this
% same image.
%
% ## Output
%
% ### Probe detection results
% A '.mat' file containing the following variables:
%
% - 'model_filename': The path to the file containing user-provided
%   measurements of the probe in the structure 'probe'. The 'probe'
%   structure is provided to this script via the output of
%   '.\CreateProbeDetectionModel.m', and so 'model_filename' is actually a
%   direct parameter of '.\CreateProbeDetectionModel.m' and an indirect
%   parameter of this script. It is retrieved from the output of
%   '.\CreateProbeDetectionModel.m' and copied to the output of this script
%   for completeness.
%
% - 'probe_detection_matches': A structure vector describing the matches
%   found from detected edges between bands of colour on the probe in the
%   image to user-provided measurements of the edges between probe colour
%   bands. Each element contains the match information for one detected
%   edge. If no suitable match was found for a given edge detected in the
%   image, the corresponding field values will be `NaN`. The fields of the
%   structure vector are as follows:
%   - index: The index of the detected edge between two bands of colour
%       on the probe. The indices of edges correspond to the order of the
%       edges along the major axis of the probe in the image, from left to
%       right.
%   - lengthAlongPCAMajorAxis: The pixel distance along the major axis of
%       the probe from the first detected edge to the current edge. For
%       more information, refer to the documentation of the `lengths`
%       output argument of `bilateralModel`.
%   - pointAbovePCAMajorAxis: A two-element row vector containing the detected
%       position of the point where the edge meets the lower border of the
%       probe in the image. (This is the point "above" the axis of the
%       probe in the sense that it has a pixel y-coordinate that is larger
%       than that of the axis at the same pixel x-coordinate.)
%   - pointBelowPCAMajorAxis: Similar to 'pointAbovePCAMajorAxis', but contains
%       the detected position of the point where the edge meets the upper
%       border of the probe in the image.
%   - matchedLengthIndex: The index of the edge between coloured bands of
%       the physical probe that is matched with the detected edge.
%       Specifically, this is the index into the user-provided measurements
%       of the probe, `probe.lengths` and `probe.widths`, where `probe` is
%       a structure input and output by '.\CreateProbeDetectionModel.m'.
%   - matchedLength: The value of `probe.lengths(matchedLengthIndex)`, where
%       `probe` is a structure input and output by
%       '.\CreateProbeDetectionModel.m'.
%   - matchedWidth: The value of `probe.widths(matchedLengthIndex)`, where
%       `probe` is a structure input and output by
%       '.\CreateProbeDetectionModel.m'.
%
% - 'probe_detection_matches_filtered': A copy of
%   'probe_detection_matches', with detected edges that were not matched to
%   probe measurements filtered out. In other words, this is a version of
%   'probe_detection_matches' without elements containing `NaN` values for
%   any fields.
%
% Additionally, the output file contains the values of all parameters in
% the first section of the script below, for reference. (Specifically,
% those listed in `parameters_list`, which should be updated if the set of
% parameters is changed.)
%
% ## References
% - M.-C. Chuang, J.-N. Hwang, K. Williams and R. Towler. "Tracking Live
%   Fish from Low-Contrast and Low-Frame-Rate Stereo Videos". IEEE
%   Transactions on Circuits and Systems for Video Technology, vol. 25, no.
%   1, pp. 167-179, Jan. 2015.
% - T. Gevers and H. Stokman. "Robust Histogram Construction from Color
%   Invariants for Object Recognition". IEEE Transactions on Pattern
%   Analysis and Machine Intelligence, vol. 26, no. 1, pp. 113-118, Jan.
%   2004.

% Bernard Llanos
% Spring 2016 research assistantship supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created August 10, 2016

%% Input data and parameters

% List of parameters to save with results
parameters_list = {
        'detection_model_filename',...
        'I_filename',...
        'rgb_sigma_filename',...
        'detectBoundingBoxesParams',...
        'detectWithinBoundingBoxesParams'
    };

% Probe detection model
detection_model_filename = 'C:\Users\llanos\Google Drive\PointProbing\Data and results\20160811_bambooSkewerProbe\20170112_bambooSkewer_orangeBlue_probeDetectionModel_bottomCamera_rect.mat';
% Image of probe in an unknown pose
I_filename = 'C:\Users\llanos\Google Drive\PointProbing\Data and results\20160811_bambooSkewerProbe\undistorted\probePrePaperOcclusion_1_b_rect.bmp';
% RGB noise parameters
rgb_sigma_filename = 'C:\Users\llanos\Google Drive\PointProbing\Data and results\20160811_bambooSkewerProbe\20160811_rgbStddev_bottomCamera.mat';

% Ask for probe's bounding region
request_bounding_region = true;

% Determination of the probe's bounding region
detectBoundingBoxesParams.noise_threshold = []; % Select automatically using Otsu's method
erosion_radius = 5;
detectBoundingBoxesParams.erosion_radius = erosion_radius;
detectBoundingBoxesParams.radius_adj = 2 * erosion_radius + 10;
detectBoundingBoxesParams.axis_distance_outlier_threshold = 3;
detectBoundingBoxesParams.dilation_radius = 4 * erosion_radius;
detectBoundingBoxesParams.region_expansion_factor_length = 1.1;
detectBoundingBoxesParams.region_expansion_factor_width = 1.5;

% Determination of refined probe colour regions
% Threshold for identifying noise pixels in final histogram backprojections
noise_threshold_final = 0.2;
% Radius for eroding images
erosion_radius_final = 2;
% Radius used to filter probe colour regions to those close to regions for
% other colours
radius_adj_final = 2 * erosion_radius_final + 4;
% Number of standard deviations from the estimate of the probe axis beyond
% which a region is determined to be distinct from the probe
axis_distance_outlier_threshold_final = 3;

% Location of probe edge points
% Search distance from the probe colour regions for band junction pixels
band_edge_distance_threshold = 2 * erosion_radius_final + 1;

% Location of probe edge endpoints
% Characteristic edge width
edge_refinement_edge_width = band_edge_distance_threshold;
% Standard deviation for the Gaussian filter applied to edge orientations
edge_refinement_angle_std = pi / 12;
% Threshold for the filtered edge image
edge_refinement_filter_threshold = 0.3;

% Matching edge endpoints to probe measurements
detected_point_alignment_outlier_threshold = 5;
color_sum_outlier_threshold = 1;
color_dominance_threshold = 1.5;
subject_gap_cost_detection = -1;
query_gap_cost_detection = 0;
% Set a low weight when the pattern of probe colours is nearly symmetrical
color_weight_detection = 0.5;

% Debugging tools
detectBoundingBoxesVerbose.display_original_image = true;
detectBoundingBoxesVerbose.display_hue_image = true;
detectBoundingBoxesVerbose.plot_global_hue_estimator = true;
detectBoundingBoxesVerbose.plot_ratio_estimators = true;
detectBoundingBoxesVerbose.display_ratio_distribution_backprojections = true;
detectBoundingBoxesVerbose.verbose_region_extraction = true;
detectBoundingBoxesVerbose.verbose_region_filtering = true;
detectBoundingBoxesVerbose.display_region_expansion = true;
plot_bounding_area_hue_estimator = false;
plot_final_ratio_estimators = false;
display_final_ratio_distribution_backprojections = false;
verbose_final_region_extraction = false;
verbose_final_region_filtering = false;
display_final_regions_colored = false;
display_band_edge_extraction = false;
verbose_edge_endpoint_extraction = false;
display_detected_model_from_image = true;
display_final_clipped_regions_colored = false;
verbose_detected_point_sequence_matching = false;
display_detected_model_matching = true;

%% Load the image containing the probe in an unknown pose

I = imread(I_filename);
image_width = size(I, 2);
image_height = size(I, 1);
image_n_channels = size(I, 3);
if image_n_channels ~= 3
    error('A 3-channel RGB image is required.')
end

%% Load the probe detection model
model_variables_required = {...
        'probe_color_distribution_resolution',...
        'probe',...
        'probe_color_distributions',...
        'probe_color_distribution_increment',...
        'model_filename'...
    };
load(detection_model_filename, model_variables_required{:});
if ~all(ismember(model_variables_required, who))
    error('One or more of the probe detection model variables is not loaded.')
end

n_colors = size(probe_color_distributions, 2);

%% Compute variables used later

% Obtain hue values
H = rgb2hue(I);

load(rgb_sigma_filename, 'rgb_sigma_polyfit');
if ~exist('rgb_sigma_polyfit', 'var')
    error('No variable called ''rgb_sigma_polyfit'' is loaded (which would contain the camera RGB noise model).')
end

I_double = im2double(I);
R = I_double(:, :, 1);
G = I_double(:, :, 2);
B = I_double(:, :, 3);

% Identify which pairs of colours can be adjacent
probe_color_pairs = [probe.colors(1:(end - 1)), probe.colors(2:end)];
probe_color_pairs = sort(probe_color_pairs, 2);
probe_color_pairs = unique(probe_color_pairs, 'rows');
n_color_pairs = size(probe_color_pairs, 1);

%% Manual bounding region selection (if enabled)
if request_bounding_region
    figure
    imshow(I);
    title('Select a bounding region')
    probe_mask_initial = roipoly;
    if isempty(probe_mask_initial)
        error('Bounding region selection aborted. Either try again, or set `request_bounding_region` to false.');
    end
else
    probe_mask_initial = detectBoundingBoxes(...
        I, probe,...
        probe_color_distribution_resolution, probe_color_distributions,...
        probe_color_distribution_increment, rgb_sigma_polyfit,...
        detectBoundingBoxesParams, detectBoundingBoxesVerbose...
    );
end

%% Compute the hue variable kernel density estimator for the bounding region

[...
    bound_color_distribution,...
    bound_color_distribution_increment...
] = hueVariableKernelDensityEstimator(...
    H, R, G, B, probe_mask_initial,...
    rgb_sigma_polyfit, probe_color_distribution_resolution...
);

if plot_bounding_area_hue_estimator
    legend_names = cell(n_colors + 1, 1); %#ok<UNRCH>
    legend_names{1} = 'Bounding area';
    for i = 1:n_colors
        legend_names{i + 1} = sprintf('Probe colour %d', i);
    end
    plotHueVariableKernelDensityEstimator(...
        bound_color_distribution_increment,...
        [bound_color_distribution, probe_color_distributions], legend_names...
    );
    title('Hue variable kernel density estimator for the initial bounding area')
end

%% Transform the image using histogram backprojection for the bounding region

% Ratio distributions with respect to the background
ratio_distributions_bounding = zeros(...
        probe_color_distribution_resolution, n_colors...
    );
for i = 1:n_colors
    ratio_distributions_bounding(:, i) = ratioDistribution(...
            probe_color_distributions(:, i), bound_color_distribution...
        );
end

if plot_final_ratio_estimators
    legend_names = cell(n_colors, 1); %#ok<UNRCH>
    for i = 1:n_colors
        legend_names{i} = sprintf('Probe colour %d', i);
    end
    plotHueVariableKernelDensityEstimator(...
        probe_color_distribution_increment, ratio_distributions_bounding, legend_names...
    );
    title('Ratio hue variable kernel density estimators for probe colours with respect to the bounding area')
end

% Histogram backprojection
ratio_distributions_backprojected_bounding = zeros(image_height, image_width, n_colors);
for i = 1:n_colors
    ratio_distributions_backprojected_bounding(:, :, i) = queryDiscretized1DFunction(...
            H, ratio_distributions_bounding(:, i), bound_color_distribution_increment...
        );
end

if display_final_ratio_distribution_backprojections
    for i = 1:n_colors %#ok<UNRCH>
        figure
        imshow(...
                ratio_distributions_backprojected_bounding(:, :, i) /...
                max(max(ratio_distributions_backprojected_bounding(:, :, i)))...
            );
        title(sprintf('Ratio distribution backprojection for probe band %d with respect to the bounding area', i))
    end
end

%% Find final regions corresponding to probe colours

[ probe_regions_final, probe_regions_bw_final] = extractBinaryRegions(...
        ratio_distributions_backprojected_bounding,...
        noise_threshold_final,...
        erosion_radius_final,...
        probe_mask_initial,...
        verbose_final_region_extraction...
    );

[
    probe_regions_final_filtered,...
    probe_regions_bw_final_filtered...
] = detectProbeBinaryRegions(...
        probe_regions_final,...
        probe_regions_bw_final,...
        probe_color_pairs,...
        radius_adj_final,...
        axis_distance_outlier_threshold_final,...
        verbose_final_region_filtering...
    );

if display_final_regions_colored
    probe_regions_bw_final_display = zeros(image_height, image_width, image_n_channels); %#ok<UNRCH>
    for i = 1:n_colors
        [ ~, peak_hue_index ] = max(probe_color_distributions(:, i));
        peak_hue = (peak_hue_index - 1) * probe_color_distribution_increment;
        peak_rgb = hsv2rgb([peak_hue, 1, 1]);
        probe_regions_bw_final_filtered_i = probe_regions_bw_final_filtered(:, :, i);
        probe_regions_bw_final_display = probe_regions_bw_final_display +...
            cat(3,...
                    peak_rgb(1) * probe_regions_bw_final_filtered_i,...
                    peak_rgb(2) * probe_regions_bw_final_filtered_i,...
                    peak_rgb(3) * probe_regions_bw_final_filtered_i...
                );
    end
    figure
    imshow(probe_regions_bw_final_display);
    title('Final detected probe regions')
end

%% Extract points along the edges between coloured bands

% Find candidate edges as places where the distance to adjacent pairs of
% colours is small
probe_color_pairs_bwdist = zeros(image_height, image_width, n_color_pairs);
for i = 1:n_color_pairs
    pair = probe_color_pairs(i, :);
    probe_regions_bw_final_filtered_1 = probe_regions_bw_final_filtered(:, :, pair(1));
    distances_color1 = bwdist(probe_regions_bw_final_filtered_1);
    distances_color1(probe_regions_bw_final_filtered_1) = Inf;
    probe_regions_bw_final_filtered_2 = probe_regions_bw_final_filtered(:, :, pair(2));
    distances_color2 = bwdist(probe_regions_bw_final_filtered_2);
    distances_color2(probe_regions_bw_final_filtered_2) = Inf;
    probe_color_pairs_bwdist(:, :, i) = max(distances_color1, distances_color2);
end
probe_color_pairs_bwdist_all = min(probe_color_pairs_bwdist, [], 3);
probe_color_pairs_bwdist_all = (probe_color_pairs_bwdist_all <= band_edge_distance_threshold);

if display_band_edge_extraction
    figure %#ok<UNRCH>
    imshow(probe_color_pairs_bwdist_all);
    title(sprintf('Euclidean distances of %g or less to adjacent probe colour regions', band_edge_distance_threshold))
end

% Refine the edges between bands
interest_points_detected = detectProbeEdgeEndpoints(...
        probe_color_pairs_bwdist_all,...
        edge_refinement_edge_width,...
        edge_refinement_angle_std,...
        edge_refinement_filter_threshold,...
        verbose_edge_endpoint_extraction...
    );

%% Organize the detected points into per-edge pairs, and filter outliers

% In contrast to probe model creation, the probe tips are not extracted,
% and therefore not expected as part of the output of `bilateralModel`. The
% `distinguish_tips` parameter is therefore set to `false`.
[...
    model_from_image_detected_pca_space,...
    image_lengths_detected,...
    model_from_image_detected_axes,...
    model_from_image_detected,...
    model_from_image_transform...
] = bilateralModel(...
    interest_points_detected, detected_point_alignment_outlier_threshold, false...
);

n_detected_band_edges = length(image_lengths_detected);

if display_detected_model_from_image
    fg = figure; %#ok<UNRCH>
    imshow(I);
    plotBilateralModel( model_from_image_detected, model_from_image_detected_axes, [image_height, image_width], [], fg );
    title('Classified detected interest points (green, red = above/below first PCA component; yellow = unmatched)');
end

if ~isempty(model_from_image_detected.unmatched)
    warning(sprintf(['Not all detected probe colour band junctions are marked with two points - One on each edge of the probe.\n',...
        'The output of probe colour band junction detection may have been quite messy.\n',...
        'Consider adjusting edge detection and edge endpoint detection parameters,\n',...
        'or increasing the outlier detection threshold used when pairing points.'])); %#ok<SPWRN>
end

%% Label the colours between detected probe band edges

% Clip the detected colour regions to the area between the probe edge
% endpoints. Extrapolate the clipping area to the tips of the probe by extending the
% line segments formed by the first and last pairs of band edges.

% Interior region
polygons_points = zeros((n_detected_band_edges + 2) * 2, 2);
polygons_points(2:(n_detected_band_edges + 1), :) = model_from_image_detected.above;
polygons_points((n_detected_band_edges + 4):(end-1), :) = flipud(model_from_image_detected.below);

% Define line segments used to extrapolate to the probe tips
line_endpoints = [
    model_from_image_detected.above(1,:), model_from_image_detected.above(2,:);
    model_from_image_detected.above(end,:), model_from_image_detected.above(end-1,:);
    model_from_image_detected.below(end,:), model_from_image_detected.below(end-1,:);
    model_from_image_detected.below(1,:), model_from_image_detected.below(2,:)
    ];

% Find intersections with the four image borders
% Top, right, bottom, left
border_coordinates = [1; image_width; image_height; 1];
border_coordinate_indices = [2; 1; 2; 1];

line_endpoints_rep = repmat(line_endpoints, 4, 1);
border_coordinates_rep = repelem(border_coordinates, 4);
border_coordinate_indices = repelem(border_coordinate_indices, 4);
line_indices = repmat((1:size(line_endpoints, 1))', 4, 1);

border_points = pointsOnLinesByCoordinates(...
    line_endpoints_rep(:, 1:2), line_endpoints_rep(:, 3:4),...
    border_coordinates_rep, border_coordinate_indices);

% Determine which intersection points to use
border_points_filter = all(isfinite(border_points), 2);
% Select intersection points in the appropriate directions
border_points_filter = border_points_filter &...
    (dot(border_points - line_endpoints_rep(:, 1:2),...
     line_endpoints_rep(:, 1:2) - line_endpoints_rep(:, 3:4), 2) > 0);
% Select intersection points within the image
border_points_filter = border_points_filter &...
    border_points(:, 1) > 0 & border_points(:, 1) <= image_width &...
    border_points(:, 2) > 0 & border_points(:, 2) <= image_height;
line_indices = line_indices(border_points_filter);
border_points = border_points(border_points_filter, :);
border_points(line_indices, :) = border_points;

% If lines crossed, replace border points by line intersections
if dot(border_points(1, :) - border_points(4, :), line_endpoints(1, 1:2) - line_endpoints(4, 1:2), 2) < 0
    intersection_point = cross(...
            cross([border_points(1, :), 1], [line_endpoints(1, 1:2), 1]),...
            cross([border_points(4, :), 1], [line_endpoints(4, 1:2), 1])...
        );
    border_points([1 4], :) = repmat(intersection_point(1:2) ./ intersection_point(3), 2, 1);
end
if dot(border_points(2, :) - border_points(3, :), line_endpoints(2, 1:2) - line_endpoints(3, 1:2), 2) < 0
    intersection_point = cross(...
            cross([border_points(2, :), 1], [line_endpoints(2, 1:2), 1]),...
            cross([border_points(3, :), 1], [line_endpoints(3, 1:2), 1])...
        );
    border_points([2 3], :) = repmat(intersection_point(1:2) ./ intersection_point(3), 2, 1);
end

polygons_points([1, n_detected_band_edges + 2, n_detected_band_edges + 3, end], :) = border_points;

% Perform the clipping
clipping_mask = roipoly(I, polygons_points(:, 1), polygons_points(:, 2));
probe_regions_bw_final_clipped = probe_regions_bw_final_filtered;

for i = 1:n_colors
    probe_regions_bw_final_clipped(:, :, i) = probe_regions_bw_final_clipped(:, :, i) & clipping_mask;
end

if display_final_clipped_regions_colored
    probe_regions_bw_final_clipped_display = zeros(image_height, image_width, image_n_channels); %#ok<UNRCH>
    for i = 1:n_colors
        [ ~, peak_hue_index ] = max(probe_color_distributions(:, i));
        peak_hue = (peak_hue_index - 1) * probe_color_distribution_increment;
        peak_rgb = hsv2rgb([peak_hue, 1, 1]);
        probe_regions_bw_final_clipped_filtered_i = probe_regions_bw_final_clipped(:, :, i);
        probe_regions_bw_final_clipped_display = probe_regions_bw_final_clipped_display +...
            cat(3,...
                    peak_rgb(1) * probe_regions_bw_final_clipped_filtered_i,...
                    peak_rgb(2) * probe_regions_bw_final_clipped_filtered_i,...
                    peak_rgb(3) * probe_regions_bw_final_clipped_filtered_i...
                );
    end
    figure
    imshow(probe_regions_bw_final_clipped_display);
    title('Final detected probe regions, clipped to detected probe edges')
end

% Find the first components of the PCA-space coordinates of the probe colour regions
probe_regions_pca_space_x = cell(n_colors, 1);
for i = 1:n_colors
    probe_regions_bw_final_clipped_i = probe_regions_bw_final_clipped(:, :, i);
    pixel_indices_i = find(probe_regions_bw_final_clipped_i);
    [y,x] = ind2sub([image_height, image_width],pixel_indices_i);
    probe_regions_pca_space_xi = [x, y, ones(length(x),1)];
    probe_regions_pca_space_xi = probe_regions_pca_space_xi / model_from_image_transform;
    % Normalization (dividing by the homogenous coordinate) is not necessary
    probe_regions_pca_space_x{i} = probe_regions_pca_space_xi(:, 1);
end

% Iterate over each probe band to find the colours dominating at either end
probe_color_scores_left = zeros(n_detected_band_edges, n_colors);
probe_color_scores_right = zeros(n_detected_band_edges, n_colors);
for i = 1:n_detected_band_edges
    if i > 1
        min_coord = mean(...
            [model_from_image_detected_pca_space.above(i-1,1),...
            model_from_image_detected_pca_space.below(i-1,1)]...
            );
    else
        min_coord = -Inf;
    end
    if i < n_detected_band_edges
        max_coord = mean(...
            [model_from_image_detected_pca_space.above(i+1,1),...
            model_from_image_detected_pca_space.below(i+1,1)]...
            );
    else
        max_coord = Inf;
    end
    center_coord = mean(...
            [model_from_image_detected_pca_space.above(i,1),...
            model_from_image_detected_pca_space.below(i,1)]...
            );
    for j = 1:n_colors
        probe_regions_pca_space_xi = probe_regions_pca_space_x{j};
        left_x = probe_regions_pca_space_xi > min_coord &...
            probe_regions_pca_space_xi <= center_coord;
        left_x = probe_regions_pca_space_xi(left_x);
        left_x = repmat(center_coord, length(left_x), 1) - left_x;
        probe_color_scores_left(i,j) = sum(left_x);
        right_x = probe_regions_pca_space_xi >= center_coord &...
            probe_regions_pca_space_xi < max_coord;
        right_x = probe_regions_pca_space_xi(right_x);
        right_x = right_x - repmat(center_coord, length(right_x), 1);
        probe_color_scores_right(i,j) = sum(right_x);
    end
end

% Pick the dominant colour, provided it dominates by a sufficient factor
for i = 1:2
    if i == 1
        scores = probe_color_scores_left;
    else
        scores = probe_color_scores_right;
    end
    [scores_max, scores_maxind] = max(scores, [], 2);
    
    % Filter out outlier small scores
    mu_scores_max = mean(scores_max);
    sigma_scores_max = std(scores_max);
    mu_scores_max_rep = repmat(mu_scores_max,n_detected_band_edges,1);
    sigma_scores_max_rep = repmat(sigma_scores_max,n_detected_band_edges,1);
    scores_max_filter = (mu_scores_max - scores_max) > color_sum_outlier_threshold * sigma_scores_max_rep;
    scores_max(scores_max_filter) = 0;

    scores_nomax = scores;
    scores_nomax(...
        sub2ind(size(scores), (1:n_detected_band_edges)', scores_maxind)...
        ) = 0;
    scores_max2 = max(scores_nomax, [], 2);
    colors = scores_max ./ scores_max2;
    colors_filter = colors >= color_dominance_threshold;
    colors(~colors_filter) = -1; % Any colour
    colors(scores_max == 0) = 0; % No colour
    colors(colors_filter) = scores_maxind(colors_filter);
    if i == 1
        probe_colors_detected_left = colors;
    else
        probe_colors_detected_right = colors;
    end
end

%% Match model extracted from the image to the user-supplied measurements of the probe
if ~exist('probe', 'var')
    error('No variable called ''probe'' is loaded (which would contain probe measurements).')
end

% Note that the probe tips should never be in the sequence of interest points
% detected in the image, because they are not adjacent to any coloured
% bands. But, in case they are, match them and then remove them from the
% matched sequence.
matching_start_index = 1;
matching_end_index = length(probe.lengths);
probe_lengths_for_matching = probe.lengths(matching_start_index:matching_end_index);
subject = [probe_lengths_for_matching, zeros(matching_end_index - matching_start_index + 1, 2)];
% Colours to the left of each band transition
subject((matching_start_index + 1):matching_end_index, 2) = probe.colors;
% Colours to the right of each band transition
subject(matching_start_index:(matching_end_index - 1), 3) = probe.colors;

query = [image_lengths_detected, probe_colors_detected_left, probe_colors_detected_right];

image_to_measured_matches_detected = matchProbeLengths(...
  subject(:, 1),...
  query(:, 1),...
  subject_gap_cost_detection,...
  query_gap_cost_detection,...
  subject(:, 2:3),...
  query(:, 2:3),...
  color_weight_detection,...
  verbose_detected_point_sequence_matching...
);

% Remove probe tips
tip_filter = image_to_measured_matches_detected == matching_start_index;
if any(tip_filter)
    warning('Removing detected interest point which was matched to the start of the probe.')
    image_to_measured_matches_detected(tip_filter) = 0;
end
tip_filter = image_to_measured_matches_detected == matching_end_index;
if any(tip_filter)
    warning('Removing detected interest point which was matched to the end of the probe.')
    image_to_measured_matches_detected(tip_filter) = 0;
end

%% Organize matching results into a structure for output

probe_widths_for_matching = probe.widths(matching_start_index:matching_end_index);
image_to_measured_matches_detected_filter = logical(image_to_measured_matches_detected);
image_to_measured_matches_detected(~image_to_measured_matches_detected_filter) = NaN;
if ~all(image_to_measured_matches_detected_filter)
    warning('Some detected interest points in the image were not matched to known probe measurements.')
    matched_lengths = nan(n_detected_band_edges, 1);
    matched_lengths(image_to_measured_matches_detected_filter) =...
        probe_lengths_for_matching(...
                image_to_measured_matches_detected(...
                        image_to_measured_matches_detected_filter...
                    )...
            );
    matched_widths = nan(n_detected_band_edges, 1);
    matched_widths(image_to_measured_matches_detected_filter) =...
        probe_widths_for_matching(...
                image_to_measured_matches_detected(...
                        image_to_measured_matches_detected_filter...
                    )...
            );
else
    matched_lengths = probe_lengths_for_matching(image_to_measured_matches_detected);
    matched_widths = probe_widths_for_matching(image_to_measured_matches_detected);
end

probe_detection_matches = struct(...
        'index', num2cell((1:n_detected_band_edges).'),...
        'lengthAlongPCAMajorAxis', num2cell(image_lengths_detected),...
        'pointAbovePCAMajorAxis', num2cell(model_from_image_detected.above, 2),...
        'pointBelowPCAMajorAxis', num2cell(model_from_image_detected.below, 2),...
        'matchedLengthIndex', num2cell(...
                image_to_measured_matches_detected + (matching_start_index - 1)...
            ),...
        'matchedLength', num2cell(matched_lengths),...
        'matchedWidth', num2cell(matched_widths)...
    );

if display_detected_model_matching
    disp('Match between detected interest points and probe measurements:') %#ok<UNRCH>
    probe_detection_matches_display =...
        rmfield(probe_detection_matches, 'lengthAlongPCAMajorAxis');
    probe_detection_matches_display =...
        rmfield(probe_detection_matches_display, 'matchedWidth');
    probe_detection_matches_display = struct2table(probe_detection_matches_display);
    disp(probe_detection_matches_display);
end

%% Save results to a file
probe_detection_matches_filtered = probe_detection_matches(image_to_measured_matches_detected_filter);
save_variables_list = [ parameters_list, {...
        'model_filename',...
        'probe_detection_matches',...
        'probe_detection_matches_filtered'
    } ];
uisave(save_variables_list,'probeDetectionResults')