%% Probe Object Detection
% Use previously estimated colour distributions and measured dimensions of
% the probe object to locate it in a test image.
%
% ## Usage
%   Modify parameters and paths to input data in the first code section
%   below, then run.
%
% ## Probe Detection Assumptions and Limitations
% - The image contains fairly saturated colours throughout, to assist with
%   hue-based colour discrimination.
% - The image only contains a single instance of the probe.
% - The probe may be occluded in one or more places. However, it can be
%   detected if at least 5 transition lines between pairs of coloured bands
%   of the probe are completely visible. The transition lines do not need
%   to be between bands that are all consecutive.
%
% ## Input
%
% ### Probe detection model
% A '.mat' file containing several variables, which is the output of
% '.\CreateProbeDetectionModel.m'. Refer to the documentation of
% '.\CreateProbeDetectionModel.m' for details.
%
% ### Image containing the probe
% An RGB image in any format that can be loaded by `imread` showing enough
% of the probe for detection to be possible. The image should have been
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
    };

% Probe detection model
detection_model_filename = 'C:\Users\Bernard\Documents\Data\20160811_bambooSkewerProbe\20160815_bambooSkewer_orangeBlue_probeDetectionModel_bottomCamera.mat';
% Image of probe in an unknown pose
I_filename = 'C:\Users\Bernard\Documents\Data\20160811_bambooSkewerProbe\original\probePrePaperOcclusion_1_b.bmp';
% RGB noise parameters
rgb_sigma_filename = 'C:\Users\Bernard\Documents\Data\20160811_bambooSkewerProbe\20160811_rgbStddev_bottomCamera.mat';

% Determination of the probe's bounding region
% Threshold for identifying noise pixels in initial histogram backprojections
noise_threshold_initial = []; % Select automatically using Otsu's method
% Radius for eroding images used to find initial bounds for the probe
erosion_radius_initial = 5;
% Radius used to filter candidate probe colour regions to those close to
% regions for other colours
radius_adj_initial = 2 * erosion_radius_initial + 10;
% Number of standard deviations from the initial estimate of the probe axis
% beyond which a region is determined to be distinct from the probe
axis_distance_outlier_threshold_initial = 3;
% Radius for dilating the candidate probe colour regions to include some of
% the background
dilation_radius_initial = 4 * erosion_radius_initial;
% Factor by which to expand oriented boxes fitted to the initial probe
% colour regions
initial_region_expansion_factor_length = 1.1;
initial_region_expansion_factor_width = 1.5;

% Determination of refined probe colour regions
% Threshold for identifying noise pixels in final histogram backprojections
noise_threshold_final = 0.2;
% Radius for eroding images
erosion_radius_final = 2;
% Radius used to filter probe colour regions to those close to regions for
% other colours
radius_adj_final = 2 * erosion_radius_final + 10;
% Number of standard deviations from the estimate of the probe axis beyond
% which a region is determined to be distinct from the probe
axis_distance_outlier_threshold_final = 3;

% Location of probe edge points
% Search distance from the probe colour regions for band junction pixels
band_edge_distance_threshold = 2 * erosion_radius_final + 1;

% Debugging tools
display_original_image = false;
display_hue_image = false;
plot_global_hue_estimator = false;
plot_initial_ratio_estimators = false;
display_initial_ratio_distribution_backprojections = false;
verbose_initial_region_extraction = false;
verbose_initial_region_filtering = false;
display_initial_region_expansion = false;
plot_bounding_area_hue_estimator = false;
plot_final_ratio_estimators = false;
display_final_ratio_distribution_backprojections = false;
verbose_final_region_extraction = false;
verbose_final_region_filtering = false;
display_final_regions_colored = true;
display_band_edge_extraction = true;

%% Load the image containing the probe in an unknown pose

I = imread(I_filename);
image_width = size(I, 2);
image_height = size(I, 1);
image_n_channels = size(I, 3);
if image_n_channels ~= 3
    error('A 3-channel RGB image is required.')
end
if display_original_image
    figure; %#ok<UNRCH>
    imshow(I);
    title('Base image');
end

%% Load the probe detection model
model_variables_required = {...
        'probe_color_distribution_resolution',...
        'probe',...
        'probe_color_distributions',...
        'probe_color_distribution_increment'...
    };
load(detection_model_filename, model_variables_required{:});
if ~all(ismember(model_variables_required, who))
    error('One or more of the probe detection model variables is not loaded.')
end

%% Compute the hue variable kernel density estimator for the image

% Obtain hue values
H = rgb2hue(I);
    
if display_hue_image
    figure %#ok<UNRCH>
    imshow(H);
    title('Hue channel of original image')
end

% Compute the hue variable kernel density estimator for the image
load(rgb_sigma_filename, 'rgb_sigma_polyfit');
if ~exist('rgb_sigma_polyfit', 'var')
    error('No variable called ''rgb_sigma_polyfit'' is loaded (which would contain the camera RGB noise model).')
end

I_double = im2double(I);
R = I_double(:, :, 1);
G = I_double(:, :, 2);
B = I_double(:, :, 3);
[...
    I_color_distribution,...
    I_color_distribution_increment...
] = hueVariableKernelDensityEstimator(...
    H, R, G, B, true(image_height, image_width),...
    rgb_sigma_polyfit, probe_color_distribution_resolution...
);

if plot_global_hue_estimator
    x = 0:I_color_distribution_increment:1; %#ok<UNRCH>
    figure
    plot(x, I_color_distribution)
    title('Hue variable kernel density estimator for the entire image')
    xlabel('Hue, \theta (range [0, 1])')
    ylabel('Density, P(\theta)')
end

%% Transform the image using histogram backprojection
n_colors = size(probe_color_distributions, 2);

% Ratio distributions with respect to the background
ratio_distributions_bg = zeros(...
        probe_color_distribution_resolution, n_colors...
    );
for i = 1:n_colors
    ratio_distributions_bg(:, i) = ratioDistribution(...
            probe_color_distributions(:, i), I_color_distribution...
        );
end

if plot_initial_ratio_estimators
    x = 0:probe_color_distribution_increment:1; %#ok<UNRCH>
    line_styles = {'-', '--', ':', '-.'};
    legend_names = cell(n_colors, 1);
    plot_colors = jet(n_colors);
    figure
    hold on
    for i = 1:n_colors
        legend_names{i} = sprintf('Probe colour %d', i);
        plot(...
                x, ratio_distributions_bg(:, i),...
                'Color', plot_colors(i, :),...
                'LineStyle', line_styles{mod(i - 1, length(line_styles)) + 1}...
            )
    end
    hold off
    legend(legend_names{:});
    title('Ratio hue variable kernel density estimators for probe colours with respect to the background')
    xlabel('Hue, \theta (range [0, 1])')
    ylabel('Density, P(\theta)')
end

% Histogram backprojection
ratio_distributions_backprojected_bg = zeros(image_height, image_width, n_colors);
for i = 1:n_colors
    ratio_distributions_backprojected_bg(:, :, i) = queryDiscretized1DFunction(...
            H, ratio_distributions_bg(:, i), I_color_distribution_increment...
        );
end

if display_initial_ratio_distribution_backprojections
    for i = 1:n_colors %#ok<UNRCH>
        figure
        imshow(...
                ratio_distributions_backprojected_bg(:, :, i) /...
                max(max(ratio_distributions_backprojected_bg(:, :, i)))...
            );
        title(sprintf('Ratio distribution backprojection for probe band %d with respect to the background', i))
    end
end

%% Find an initial bounding area for the probe

[ probe_regions_initial, probe_regions_bw_initial] = extractBinaryRegions(...
        ratio_distributions_backprojected_bg,...
        noise_threshold_initial,...
        erosion_radius_initial,...
        [],...
        verbose_initial_region_extraction...
    );

[
    ~,...
    probe_regions_bw_initial_filtered...
] = detectProbeBinaryRegions(...
        probe_regions_initial,...
        probe_regions_bw_initial,...
        radius_adj_initial,...
        axis_distance_outlier_threshold_initial,...
        verbose_initial_region_filtering...
    );

% Obtain an upper bound on the probe's area
probe_regions_bw_initial_all = any(probe_regions_bw_initial_filtered, 3);
if display_initial_region_expansion
    figure %#ok<UNRCH>
    imshow(probe_regions_bw_initial_all);
    title('Initial probe colour regions prior to expansion')
end
disk = strel('disk', dilation_radius_initial);
probe_regions_bw_initial_all = imdilate(probe_regions_bw_initial_all, disk);
% Find oriented bounding boxes
initial_region_bounds = regionprops(probe_regions_bw_initial_all, {...
    'Centroid', 'Orientation', 'MajorAxisLength', 'MinorAxisLength'...
    });
n_initial_regions = length(initial_region_bounds);
initial_region_rect_corners = cell(n_initial_regions, 1);
for i = 1:n_initial_regions
    half_length = initial_region_bounds(i).MajorAxisLength * initial_region_expansion_factor_length / 2;
    half_width = initial_region_bounds(i).MinorAxisLength * initial_region_expansion_factor_width / 2;
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
probe_mask_initial = false(image_height, image_width);
for i = 1:n_initial_regions
    probe_mask_initial = probe_mask_initial | roipoly(...
            I, initial_region_rect_corners{i}(:, 1), initial_region_rect_corners{i}(:, 2)...
        );
end
if display_initial_region_expansion
    figure %#ok<UNRCH>
    probe_regions_bw_initial_all_display = repmat(probe_regions_bw_initial_all, 1, 1, 3);
    probe_regions_bw_initial_all_display(:, :, 1) =...
        probe_regions_bw_initial_all_display(:, :, 1) |...
        imdilate(bwperim(probe_mask_initial), strel('disk', 2));
    imshow(double(probe_regions_bw_initial_all_display));
    title('Initial probe colour regions after dilation, with outline of initial mask')
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
    x = 0:bound_color_distribution_increment:1; %#ok<UNRCH>
    figure
    plot(x, bound_color_distribution)
    title('Hue variable kernel density estimator for the initial bounding area')
    xlabel('Hue, \theta (range [0, 1])')
    ylabel('Density, P(\theta)')
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
    x = 0:probe_color_distribution_increment:1; %#ok<UNRCH>
    line_styles = {'-', '--', ':', '-.'};
    legend_names = cell(n_colors, 1);
    plot_colors = jet(n_colors);
    figure
    hold on
    for i = 1:n_colors
        legend_names{i} = sprintf('Probe colour %d', i);
        plot(...
                x, ratio_distributions_bounding(:, i),...
                'Color', plot_colors(i, :),...
                'LineStyle', line_styles{mod(i - 1, length(line_styles)) + 1}...
            )
    end
    hold off
    legend(legend_names{:});
    title('Ratio hue variable kernel density estimators for probe colours with respect to the bounding area')
    xlabel('Hue, \theta (range [0, 1])')
    ylabel('Density, P(\theta)')
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

% Identify which pairs of colours can be adjacent
probe_color_pairs = [probe.colors(1:(end - 1)), probe.colors(2:end)];
probe_color_pairs = sort(probe_color_pairs, 2);
probe_color_pairs = unique(probe_color_pairs, 'rows');
n_color_pairs = size(probe_color_pairs, 1);

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