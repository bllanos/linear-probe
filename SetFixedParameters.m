%% Set Fixed Parameters
% Set values of parameters for probe detection and localization that
% seldomly need to be changed.
%
% ## Usage
%   Modify the parameters in the first two code sections below, as desired.
%   This script exists just to deduplicate code, and will be called by
%   other scripts.

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created January 18, 2018

%% Probe detection parameters

% Ask for probe's bounding region
request_bounding_region = false;

% Load variable kernel density estimators (true), or Gaussian density
% estimators (false) for the probe colours
use_kernel_estimators = true;

% Determination of the probe's bounding region
uniform_background_initial = true;
erosion_radius_initial = 5;
detectBoundingBoxesParams.saturation_threshold = 0.3;
detectBoundingBoxesParams.erosion_radius = erosion_radius_initial;
detectBoundingBoxesParams.radius_adj = 2 * erosion_radius_initial + 4;
detectBoundingBoxesParams.axis_distance_outlier_threshold = 3;
detectBoundingBoxesParams.dilation_radius = 4 * erosion_radius_initial;
detectBoundingBoxesParams.region_expansion_factor_length = 1.1;
detectBoundingBoxesParams.region_expansion_factor_width = 1.5;

% Determination of refined probe colour regions
uniform_background_final = true;
erosion_radius_final = 2;
detectWithinBoundingBoxParams.saturation_threshold = 0.15;
detectWithinBoundingBoxParams.erosion_radius = erosion_radius_final;
detectWithinBoundingBoxParams.radius_adj = 2 * erosion_radius_final + 4;
detectWithinBoundingBoxParams.axis_distance_outlier_threshold = 3;

% Location of probe edge points
band_edge_distance_threshold = 2 * erosion_radius_final + 1;
detectWithinBoundingBoxParams.band_edge_distance_threshold = band_edge_distance_threshold;

% Location of probe edge endpoints
detectWithinBoundingBoxParams.edge_refinement_edge_width = band_edge_distance_threshold;
detectWithinBoundingBoxParams.edge_refinement_angle_std = pi / 12;
detectWithinBoundingBoxParams.edge_refinement_filter_threshold = 0.3;

% Matching edge endpoints to general probe structure
detected_point_alignment_outlier_threshold = 5;

% Matching edge endpoints to probe measurements
matching.subject_gap_cost_detection = 0;
matching.query_gap_cost_detection = 0;
matching.direction_threshold = 1.25;
matching.alignment_inlier_threshold = 0.75;

% Debugging tools
detectBoundingBoxesVerbose.display_original_image = false;
detectBoundingBoxesVerbose.display_hue_image = true;
detectBoundingBoxesVerbose.display_saturation_image = true;
detectBoundingBoxesVerbose.plot_hue_estimator = true;
detectBoundingBoxesVerbose.plot_hue_classifier = true;
detectBoundingBoxesVerbose.display_distribution_backprojections = false;
detectBoundingBoxesVerbose.display_binary_images = false;
detectBoundingBoxesVerbose.verbose_region_filtering = false;
detectBoundingBoxesVerbose.display_region_expansion = true;

detectWithinBoundingBoxVerbose.display_hue_image = false;
detectWithinBoundingBoxVerbose.display_saturation_image = true;
detectWithinBoundingBoxVerbose.plot_hue_estimator = true;
detectWithinBoundingBoxVerbose.plot_hue_classifier = true;
detectWithinBoundingBoxVerbose.display_distribution_backprojections = false;
detectWithinBoundingBoxVerbose.display_binary_images = false;
detectWithinBoundingBoxVerbose.verbose_region_filtering = false;
detectWithinBoundingBoxVerbose.display_regions_colored = false;
detectWithinBoundingBoxVerbose.display_band_edge_extraction = false;
detectWithinBoundingBoxVerbose.verbose_edge_endpoint_extraction = true;

verbose.display_detected_model_from_image = true;
verbose.display_final_clipped_regions_colored = true;
verbose.verbose_detected_point_sequence_matching = false;
verbose.display_detected_model_matching = true;
verbose.warnings = true;

%% Probe localization parameters

% Linear Probe Localization
% Error convergence threshold for linear probe estimation
linear_convergence_threshold = 0.01;
% Normalize lengths when estimating a 1D homography between the probe and
% its image
normalize_homography1D = true;

% Enable or disable nonlinear iterative refinement
enable_nonlinear_estimation = true;

% Debugging tools
verbose.verbose_linear_estimation = false; % Requires `I_filename` to be valid
verbose.display_linear_estimation = false; % Requires `I_filename` to be valid
verbose.verbose_nonlinear_estimation = false;
verbose.display_nonlinear_estimation = true; % Requires `I_filename` to be valid
verbose.display_axis_points = true; % Requires `I_filename` to be valid

%% Pack probe detection parameters

detectionParams.request_bounding_region = request_bounding_region;
detectionParams.uniform_background_initial = uniform_background_initial;
detectionParams.detectBoundingBoxesParams = detectBoundingBoxesParams;
detectionParams.uniform_background_final = uniform_background_final;
detectionParams.detectWithinBoundingBoxParams = detectWithinBoundingBoxParams;
detectionParams.detected_point_alignment_outlier_threshold = detected_point_alignment_outlier_threshold;
detectionParams.matching = matching;

verbose.detectBoundingBoxesVerbose = detectBoundingBoxesVerbose;
verbose.detectWithinBoundingBoxVerbose = detectWithinBoundingBoxVerbose;

%% Pack probe localization parameters

localizationParams.linear_convergence_threshold = linear_convergence_threshold;
localizationParams.normalize_homography1D = normalize_homography1D;
localizationParams.enable_nonlinear_estimation = enable_nonlinear_estimation;