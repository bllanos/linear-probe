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
        'request_bounding_region',...
        'detectBoundingBoxesParams',...
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

% Determination of the probe's bounding region
erosion_radius_initial = 5;
detectBoundingBoxesParams.erosion_radius = erosion_radius_initial;
detectBoundingBoxesParams.radius_adj = 2 * erosion_radius_initial + 10;
detectBoundingBoxesParams.axis_distance_outlier_threshold = 3;
detectBoundingBoxesParams.dilation_radius = 4 * erosion_radius_initial;
detectBoundingBoxesParams.region_expansion_factor_length = 1.1;
detectBoundingBoxesParams.region_expansion_factor_width = 1.5;

% Determination of refined probe colour regions
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
subject_gap_cost_detection = 0;
query_gap_cost_detection = 0;
direction_threshold = 1.25;
alignment_inlier_threshold = 0.75;

% Debugging tools
detectBoundingBoxesVerbose.display_original_image = false;
detectBoundingBoxesVerbose.display_hue_image = false;
detectBoundingBoxesVerbose.plot_hue_estimator = true;
detectBoundingBoxesVerbose.plot_hue_classifier = true;
detectBoundingBoxesVerbose.display_distribution_backprojections = true;
detectBoundingBoxesVerbose.display_binary_images = true;
detectBoundingBoxesVerbose.verbose_region_filtering = false;
detectBoundingBoxesVerbose.display_region_expansion = true;

detectWithinBoundingBoxVerbose.display_hue_image = false;
detectWithinBoundingBoxVerbose.plot_hue_estimator = true;
detectWithinBoundingBoxVerbose.plot_hue_classifier = true;
detectWithinBoundingBoxVerbose.display_distribution_backprojections = true;
detectWithinBoundingBoxVerbose.display_binary_images = true;
detectWithinBoundingBoxVerbose.verbose_region_filtering = false;
detectWithinBoundingBoxVerbose.display_regions_colored = false;
detectWithinBoundingBoxVerbose.display_band_edge_extraction = false;
detectWithinBoundingBoxVerbose.verbose_edge_endpoint_extraction = false;

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
        'probe',...
        'probe_color_distributions_kernel',...
        'model_filename'...
    };
load(detection_model_filename, model_variables_required{:});
if ~all(ismember(model_variables_required, who))
    error('One or more of the probe detection model variables is not loaded.')
end

load(rgb_sigma_filename, 'rgb_sigma_polyfit');
if ~exist('rgb_sigma_polyfit', 'var')
    error('No variable called ''rgb_sigma_polyfit'' is loaded (which would contain the camera RGB noise model).')
end

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
        probe_color_distributions_kernel, rgb_sigma_polyfit,...
        detectBoundingBoxesParams, detectBoundingBoxesVerbose...
    );
end

%% Detect interest points and refined probe colour regions within the bounding box

[...
    interest_points_detected, probe_regions_bw_final_filtered...
] = detectWithinBoundingBox(...
    I, probe, probe_mask_initial,...
    probe_color_distributions_kernel, rgb_sigma_polyfit,...
    detectWithinBoundingBoxParams, detectWithinBoundingBoxVerbose...
    );

if isempty(interest_points_detected)
    error('No interest points detected.')
end

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

[...
    probe_colors_detected_left, probe_colors_detected_right...
] = labelDetectedColors(...
    I, probe_color_distributions_kernel,...
    probe_regions_bw_final_filtered,...
    model_from_image_detected_pca_space, model_from_image_detected,...
    model_from_image_transform, model_from_image_detected_axes(1, :),...
    display_final_clipped_regions_colored...
);


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

image_to_measured_matches_detected = matchProbeLengthsRandom(...
  subject(:, 1),...
  query(:, 1),...
  subject_gap_cost_detection,...
  query_gap_cost_detection,...
  subject(:, 2:3),...
  query(:, 2:3),...
  direction_threshold,...
  alignment_inlier_threshold,...
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