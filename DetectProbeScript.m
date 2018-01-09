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
%   located in 3D space if at least 3 edges between pairs of coloured bands
%   of the probe are completely visible. (Localization is not handled by
%   this script, but by 'LocalizeProbeScript.m'.) The edges do not need to
%   be between bands that are all consecutive.
%
% ## Input
%
% ### Probe detection model
% A '.mat' file containing several variables, which is the output of
% 'CreateProbeDetectionModel.m'. Refer to the documentation of
% 'CreateProbeDetectionModel.m' for details.
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
% Required if at least one of the parameters 'uniform_background_initial',
% or 'uniform_background_final' is true.
%
% A '.mat' file containing a 'rgb_sigma_polyfit' variable, as output by the
% script '.\EstimateRGBStandardDeviations.m'. 'rgb_sigma_polyfit' describes
% the variation in RGB channel standard deviations with RGB values in the
% image. This information should be computed from images taken
% under the same conditions and with the same camera parameters as the
% image in which the probe is to be detected, if not computed from this
% same image.
%
% 'rgb_sigma_polyfit' is used to compute a colour distribution for the
% image background. Instead, if 'uniform_background_initial' is true, a
% uniform distribution is used during the first pass. Likewise, if
% 'uniform_background_final' is true, a uniform distribution is used during
% the second pass.
%
% ## Output
%
% ### Probe detection results
% A '.mat' file containing the following variables:
%
% - 'model_filename': The path to the file containing user-provided
%   measurements of the probe in the structure 'probe'. The 'probe'
%   structure is provided to this script via the output of
%   'CreateProbeDetectionModel.m', and so 'model_filename' is actually a
%   direct parameter of 'CreateProbeDetectionModel.m' and an indirect
%   parameter of this script. It is retrieved from the output of
%   'CreateProbeDetectionModel.m', and copied to the output of this script
%   for completeness.
%
% - 'probe_detection_matches': The `matches` output argument of
%   'detectProbe()'. Refer to the documentation of 'detectProbe.m' for
%   details.
%
% - 'probe_detection_matches_filtered': A copy of
%   'probe_detection_matches', with detected edges that were not matched to
%   probe measurements filtered out. The `matches_filtered` output argument
%   of 'detectProbe()'. Refer to the documentation of 'detectProbe.m' for
%   details.
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
        'request_bounding_region',...
        'use_kernel_estimators',...
        'uniform_background_initial',...
        'detectBoundingBoxesParams',...
        'uniform_background_final',...
        'detectWithinBoundingBoxParams',...
        'detected_point_alignment_outlier_threshold',...
        'subject_gap_cost_detection',...
        'query_gap_cost_detection',...
        'direction_threshold',...
        'alignment_inlier_threshold'
    };

% Probe detection model
detection_model_filename = 'C:\Users\llanos\Google Drive\PointProbing\Data and results\20170605_bambooSkewerProbe_colorClassificationChanges\20170605_probeDetectionModel_2_normalization.mat';
% Image of probe in an unknown pose
I_filename = 'C:\Users\llanos\Google Drive\PointProbing\Data and results\20160811_bambooSkewerProbe\20170403_CMPUT615_Demo_data\Data\probeUnderLotus_1_b_rect.bmp';
% RGB noise parameters
rgb_sigma_filename = 'C:\Users\llanos\Google Drive\PointProbing\Data and results\20160811_bambooSkewerProbe\20170403_CMPUT615_Demo_data\Output\1_rgbStddev.mat';

% Ask for probe's bounding region
request_bounding_region = false;

% Load variable kernel density estimators (true), or Gaussian density
% estimators (false) for the probe colours
use_kernel_estimators = true;

% Determination of the probe's bounding region
uniform_background_initial = false;
erosion_radius_initial = 5;
detectBoundingBoxesParams.erosion_radius = erosion_radius_initial;
detectBoundingBoxesParams.radius_adj = 2 * erosion_radius_initial + 10;
detectBoundingBoxesParams.axis_distance_outlier_threshold = 3;
detectBoundingBoxesParams.dilation_radius = 4 * erosion_radius_initial;
detectBoundingBoxesParams.region_expansion_factor_length = 1.1;
detectBoundingBoxesParams.region_expansion_factor_width = 1.5;

% Determination of refined probe colour regions
uniform_background_final = false;
erosion_radius_final = 2;
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
detectBoundingBoxesVerbose.display_hue_image = false;
detectBoundingBoxesVerbose.plot_hue_estimator = true;
detectBoundingBoxesVerbose.plot_hue_classifier = true;
detectBoundingBoxesVerbose.display_distribution_backprojections = false;
detectBoundingBoxesVerbose.display_binary_images = false;
detectBoundingBoxesVerbose.verbose_region_filtering = false;
detectBoundingBoxesVerbose.display_region_expansion = true;

detectWithinBoundingBoxVerbose.display_hue_image = false;
detectWithinBoundingBoxVerbose.plot_hue_estimator = true;
detectWithinBoundingBoxVerbose.plot_hue_classifier = true;
detectWithinBoundingBoxVerbose.display_distribution_backprojections = false;
detectWithinBoundingBoxVerbose.display_binary_images = false;
detectWithinBoundingBoxVerbose.verbose_region_filtering = false;
detectWithinBoundingBoxVerbose.display_regions_colored = false;
detectWithinBoundingBoxVerbose.display_band_edge_extraction = true;
detectWithinBoundingBoxVerbose.verbose_edge_endpoint_extraction = true;

verbose.display_detected_model_from_image = true;
verbose.display_final_clipped_regions_colored = true;
verbose.verbose_detected_point_sequence_matching = false;
verbose.display_detected_model_matching = true;
verbose.warnings = true;

%% Load the image containing the probe in an unknown pose
I = imread(I_filename);

%% Load the probe detection model
model_variables_required = { 'probe', 'model_filename' };
if use_kernel_estimators
    model_variables_required(end + 1) = {'probe_color_distributions_kernel'};
else
    model_variables_required(end + 1) = {'probe_color_distributions_gaussian'};
end
load(detection_model_filename, model_variables_required{:});
if ~all(ismember(model_variables_required, who))
    error('One or more of the probe detection model variables is not loaded.')
end
if use_kernel_estimators
    probe_color_distributions = probe_color_distributions_kernel;
else
    probe_color_distributions = probe_color_distributions_gaussian;
end

if ~uniform_background_initial || ~uniform_background_final
    load(rgb_sigma_filename, 'rgb_sigma_polyfit');
    if ~exist('rgb_sigma_polyfit', 'var')
        error('No variable called ''rgb_sigma_polyfit'' is loaded (which would contain the camera RGB noise model).')
    end
else
    rgb_sigma_polyfit = [];
end

%% Pack parameters

params.request_bounding_region = request_bounding_region;
params.uniform_background_initial = uniform_background_initial;
params.detectBoundingBoxesParams = detectBoundingBoxesParams;
params.uniform_background_final = uniform_background_final;
params.detectWithinBoundingBoxParams = detectWithinBoundingBoxParams;
params.detected_point_alignment_outlier_threshold = detected_point_alignment_outlier_threshold;
params.matching = matching;

verbose.detectBoundingBoxesVerbose = detectBoundingBoxesVerbose;
verbose.detectWithinBoundingBoxVerbose = detectWithinBoundingBoxVerbose;

%% Detection

[ probe_detection_matches_filtered, probe_detection_matches ] = detectProbe(...
       I, probe, probe_color_distributions, rgb_sigma_polyfit, params, verbose...
);

%% Save results to a file
save_variables_list = [ parameters_list, {...
        'model_filename',...
        'probe_detection_matches',...
        'probe_detection_matches_filtered'
    } ];
uisave(save_variables_list,'probeDetectionResults');