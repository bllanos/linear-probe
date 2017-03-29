function [...
    interest_points, probe_regions_bw_filtered...
    ] = detectWithinBoundingBox(...
        I, probe, mask, probe_color_distribution_resolution,...
        probe_color_distributions,...
        probe_color_distribution_increment, rgb_sigma_polyfit,...
        params, varargin...
    )
%DETECTWITHINBOUNDINGBOX Find the probe within a bounding region
%
% ## Syntax
% [ interest_points ] = detectWithinBoundingBox(...
%       I, probe, mask, probe_color_distribution_resolution,...
%       probe_color_distributions,...
%       probe_color_distribution_increment, rgb_sigma_polyfit,...
%       params [, verbose]...
% )
% [...
%     interest_points, probe_regions_bw_filtered...
% ] = detectWithinBoundingBox(...
%       I, probe, mask, probe_color_distribution_resolution,...
%       probe_color_distributions,...
%       probe_color_distribution_increment, rgb_sigma_polyfit,...
%       params [, verbose]...
% )
%
% ## Description
% [ interest_points ] = detectWithinBoundingBox(...
%       I, probe, mask, probe_color_distribution_resolution,...
%       probe_color_distributions,...
%       probe_color_distribution_increment, rgb_sigma_polyfit,...
%       params [, verbose]...
% )
%   Returns candidate probe colour band junction endpoints within the
%   bounding region.
%
% [...
%     interest_points, probe_regions_bw_filtered...
% ] = detectWithinBoundingBox(...
%       I, probe, mask, probe_color_distribution_resolution,...
%       probe_color_distributions,...
%       probe_color_distribution_increment, rgb_sigma_polyfit,...
%       params [, verbose]...
% )
%   Additionally returns the binary images defining the probe colour bands
%   detected within the bounding region.
%
% ## Input Arguments
%
% I -- Image containing the probe
%   An RGB image containing the probe.
%
% probe -- Probe measurements
%   Refer to the documentation of './CreateProbeDetectionModel.m' for
%   details.
%
% mask -- Probe bounding region
%   A binary image with the same dimensions as `I` (but with only one
%   channel, not three), where the nonzero pixels ideally define regions
%   containing the colour bands of the probe.
%
% probe_color_distribution_resolution -- Probe colour estimator sample count
%   The number of equally-spaced samples in the range of hue values from 0
%   (inclusive) to 1 (inclusive) at which the variable kernel density
%   estimators for probe colour bands have been evaluated.
%
% probe_color_distributions -- Probe colour estimators
%   Discretized variable kernel density estimators of image hue values
%   corresponding to the different coloured bands on the probe, in the same
%   order (starting from the active tip of the probe). The i-th column of
%   this 2D array stores the estimator for the i-th colour class of probe
%   segments.
%
% probe_color_distribution_increment -- Probe colour estimator sample spacing
%   A scalar equal to the spacing between the samples of hue values in the
%   range [0, 1] at which the variable kernel density estimators have been
%   evaluated to produce 'probe_band_color_distributions'.
%
% rgb_sigma_polyfit -- Camera RGB noise model
%   An array describing the variation in RGB channel standard deviations
%   with RGB values in the image. This information should be computed from
%   images taken under the same conditions and with the same camera
%   parameters as the image (`I`) in which the probe is to be detected, if
%   not computed from this same image.
%
%   Refer to the documentation of './EstimateRGBStandardDeviations.m' for
%   details.
%
% params -- Fixed parameters
%   Parameters that should be stable across a variety of input images and
%   probe models. `params` is a structure containing the following fields:
%   - noise_threshold: Threshold for identifying noise pixels in
%       histogram backprojections. If empty (`[]`), a threshold will be
%       selected automatically using Otsu's method.
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
%   - band_edge_distance_threshold: Search distance from detected probe
%       colour regions for pixels on the edges between adjacent colour
%       regions.
%   - edge_refinement_edge_width: Characteristic ("typical") width of the
%       edges in the image between adjacent colour regions on the probe.
%   - edge_refinement_angle_std: Standard deviation parameter for the
%       Gaussian model of the orientations of edges between adjacent colour
%       regions on the probe. The image is filtered, by
%       'detectProbeEdgeEndpoints()' using a kernel constructed with this
%       parameter, and with the `edge_refinement_edge_width` parameter.
%   - edge_refinement_filter_threshold: The threshold applied to the
%       filtered image within 'detectProbeEdgeEndpoints()'.
%
% verbose -- Debugging flags
%   If recognized fields of `verbose` are true, corresponding graphical
%   output will be generated for debugging purposes.
%
%   All debugging flags default to false if `verbose` is not passed.
%
% ## Output Arguments
%
% interest_points -- Detected probe edge endpoints
%   A two column array containing the image x and y-coordinates (in the
%   first and second columns, respectively) of candidate endpoints of the
%   edges between the coloured bands of the probe.
%
%   Ideally, each edge will be represented by two endpoints in this array.
%   In practice, points will be detected along the probe's contour, far
%   from the actual edges between coloured bands, and towards the interior
%   of the probe, rather than where the edges meet the probe's contour.
%
% probe_regions_bw_filtered -- Detected probe colour bands
%   A three-dimensional logical array of size image_height x image_width x
%   n_colors. `probe_regions_bw_filtered(:, :, i)` stores the binary image
%   representing where the i-th probe colour was detected within the area
%   defined by `mask`. (`probe_regions_bw_filtered(:, :, i)` may therefore
%   define several disconnected regions, if there are multiple patches of
%   the i-th colour in the bounding mask.)
%
% See also hueVariableKernelDensityEstimator, extractBinaryRegions, detectProbeBinaryRegions, detectProbeEdgeEndpoints

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created March 28, 2017

%% Parse input arguments

nargoutchk(1,2);
narginchk(8,9);

if ~isempty(varargin)
    verbose = varargin{1};
    plot_bounding_area_hue_estimator = verbose.plot_bounding_area_hue_estimator;
    plot_ratio_estimators = verbose.plot_ratio_estimators;
    display_ratio_distribution_backprojections = verbose.display_ratio_distribution_backprojections;
    verbose_region_extraction = verbose.verbose_region_extraction;
    verbose_region_filtering = verbose.verbose_region_filtering;
    display_regions_colored = verbose.display_regions_colored;
    display_band_edge_extraction = verbose.display_band_edge_extraction;
    verbose_edge_endpoint_extraction = verbose.verbose_edge_endpoint_extraction;
else
    plot_bounding_area_hue_estimator = false;
    plot_ratio_estimators = false;
    display_ratio_distribution_backprojections = false;
    verbose_region_extraction = false;
    verbose_region_filtering = false;
    display_regions_colored = false;
    display_band_edge_extraction = false;
    verbose_edge_endpoint_extraction = false;
end

image_width = size(I, 2);
image_height = size(I, 1);
image_n_channels = size(I, 3);
if image_n_channels ~= 3
    error('A 3-channel RGB image is required.')
end

%% Compute the hue variable kernel density estimator for the bounding region

% Obtain hue values
H = rgb2hue(I);

I_double = im2double(I);
R = I_double(:, :, 1);
G = I_double(:, :, 2);
B = I_double(:, :, 3);

[...
    bound_color_distribution,...
    bound_color_distribution_increment...
] = hueVariableKernelDensityEstimator(...
    H, R, G, B, mask,...
    rgb_sigma_polyfit, probe_color_distribution_resolution...
);

n_colors = size(probe_color_distributions, 2);

if plot_bounding_area_hue_estimator
    legend_names = cell(n_colors + 1, 1);
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

if plot_ratio_estimators
    legend_names = cell(n_colors, 1);
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

if display_ratio_distribution_backprojections
    for i = 1:n_colors
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
        params.noise_threshold,...
        params.erosion_radius,...
        mask,...
        verbose_region_extraction...
    );

% Identify which pairs of colours can be adjacent
probe_color_pairs = [probe.colors(1:(end - 1)), probe.colors(2:end)];
probe_color_pairs = sort(probe_color_pairs, 2);
probe_color_pairs = unique(probe_color_pairs, 'rows');
n_color_pairs = size(probe_color_pairs, 1);

[
    ~,...
    probe_regions_bw_filtered...
] = detectProbeBinaryRegions(...
        probe_regions_final,...
        probe_regions_bw_final,...
        probe_color_pairs,...
        params.radius_adj,...
        params.axis_distance_outlier_threshold,...
        verbose_region_filtering...
    );

if display_regions_colored
    probe_regions_bw_final_display = zeros(image_height, image_width, image_n_channels);
    for i = 1:n_colors
        [ ~, peak_hue_index ] = max(probe_color_distributions(:, i));
        peak_hue = (peak_hue_index - 1) * probe_color_distribution_increment;
        peak_rgb = hsv2rgb([peak_hue, 1, 1]);
        probe_regions_bw_final_filtered_i = probe_regions_bw_filtered(:, :, i);
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
    probe_regions_bw_final_filtered_1 = probe_regions_bw_filtered(:, :, pair(1));
    distances_color1 = bwdist(probe_regions_bw_final_filtered_1);
    distances_color1(probe_regions_bw_final_filtered_1) = Inf;
    probe_regions_bw_final_filtered_2 = probe_regions_bw_filtered(:, :, pair(2));
    distances_color2 = bwdist(probe_regions_bw_final_filtered_2);
    distances_color2(probe_regions_bw_final_filtered_2) = Inf;
    probe_color_pairs_bwdist(:, :, i) = max(distances_color1, distances_color2);
end
probe_color_pairs_bwdist_all = min(probe_color_pairs_bwdist, [], 3);
probe_color_pairs_bwdist_all = (probe_color_pairs_bwdist_all <= params.band_edge_distance_threshold);

if display_band_edge_extraction
    figure
    imshow(probe_color_pairs_bwdist_all);
    title(sprintf('Euclidean distances of %g or less to adjacent probe colour regions', params.band_edge_distance_threshold))
end

% Refine the edges between bands
interest_points = detectProbeEdgeEndpoints(...
        probe_color_pairs_bwdist_all,...
        params.edge_refinement_edge_width,...
        params.edge_refinement_angle_std,...
        params.edge_refinement_filter_threshold,...
        verbose_edge_endpoint_extraction...
    );

end