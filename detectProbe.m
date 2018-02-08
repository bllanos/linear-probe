function [ matches_filtered, matches ] = detectProbe(...
        I, probe, probe_color_distributions, rgb_sigma_polyfit, params, varargin...
    )
% DETECTPROBE  Locate the probe object in an image
%
% ## Syntax
% matches_filtered = detectProbe(...
%       I, probe, probe_color_distributions, rgb_sigma_polyfit, params [, verbose]...
% )
% [ matches_filtered, matches ] = detectProbe(...
%       I, probe, probe_color_distributions, rgb_sigma_polyfit, params [, verbose]...
% )
%
% ## Description
% matches_filtered = detectProbe( I, probe, probe_color_distributions, params [, verbose] )
%   Returns structure vectors describing the probe's location in the image.
% [ matches_filtered, matches ] = detectProbe(...
%       I, probe, probe_color_distributions, params [, verbose]...
% )
%   Additionally returns a structure array describing all points detected
%   in the image, not only those which were identified with the probe.
%
% ## Input Arguments
%
% I -- Image containing the probe
%   An RGB image containing the probe.
%
% probe -- Probe measurements
%   Refer to the documentation of 'CreateProbeDetectionModel.m' for
%   details.
%
% probe_color_distributions -- Probe colour estimators
%   Discretized density estimators of image hue values corresponding to the
%   different coloured bands on the probe, in the same order (starting from
%   the active tip of the probe). The i-th column of this 2D array stores
%   the estimator for the i-th colour class of probe segments.
%
% rgb_sigma_polyfit -- Camera RGB noise model
%   A parameter passed indirectly to 'extractBinaryRegions()'. As described
%   in the documentation of 'extractBinaryRegions()'. An empty array can be
%   passed instead, if `params.uniform_background_initial` and
%   `params.uniform_background_final` are both `false`.
%
% params -- Fixed parameters
%   Parameters that should be stable across a variety of input images and
%   probe models. `params` is a structure containing the following fields:
%   - request_bounding_region: A logical scalar indicating whether to open
%     a figure and have the user draw a bounding polygon for the probe
%     (`true`), or to automatically detect the region(s) containing the
%     probe (`false`).
%   - uniform_background_initial: If `true`, assume a uniform colour
%     distribution for the image background when finding the bounding
%     region(s) containing the probe. `uniform_background_initial`
%     determines whether the `rgb_sigma_polyfit` parameter of
%     'extractBinaryRegions()' is passed empty. Refer to the documentation
%     of 'extractBinaryRegions.m' for details.
%   - detectBoundingBoxesParams: A structure of fixed parameters passed to
%     'detectBoundingBoxes()'. Refer the documentation of
%     'detectBoundingBoxes.m' for details. This field need not exist if
%     `request_bounding_region` is `true`.
%   - uniform_background_final: Similar to `uniform_background_initial`,
%     but influences the detection of the probe within its bounding
%     regions.
%   - detectWithinBoundingBoxParams: A structure of fixed parameters passed
%      to 'detectWithinBoundingBox()'. Refer the documentation of
%     'detectWithinBoundingBox.m' for details.
%   - detected_point_alignment_outlier_threshold: The
%     `point_alignment_outlier_threshold` argument of 'bilateralModel()',
%     which controls the pairing of points detected above and below the
%     probe's line of symmetry.
%   - matching: A structure containing parameters for data association
%     between the measurements of the probe, and the edges between coloured
%     bands of the probe detected in the image. `matching` has the
%     following fields, all of which are used as parameters of
%     'matchProbeLengthsRandom()'. Refer to the documentation of
%     'matchProbeLengthsRandom.m' for details.
%     - subject_gap_cost_detection: The `subject_gap_cost` parameter.
%     - query_gap_cost_detection: The `query_gap_cost` parameter.
%     - direction_threshold: The `thresholds.direction_threshold` parameter.
%     - inlier_threshold: The `thresholds.inlier_threshold`
%       parameter.
%     - local_alignment_threshold: The
%       `thresholds.local_alignment_threshold` parameter.
%
% verbose -- Debugging flags
%   If recognized fields of `verbose` are true, corresponding graphical or
%   console output will be generated for debugging purposes. Note that
%   there are two fields which are structures of debugging flags:
%   `detectBoundingBoxesVerbose`, and `detectWithinBoundingBoxVerbose`,
%   which become the `verbose` arguments of 'detectBoundingBoxes()' and
%   'detectWithinBoundingBox()', respectively. Both of these fields are
%   optional.
%
%   All debugging flags default to false if `verbose` is not passed.
%
% ## Output Arguments
%
% matches_filtered -- Detected and validated edges between coloured bands
%   A cell vector of structure vectors. The i-th cell contains a copy of
%   the i-th column of `matches`, with detected edges that were not matched
%   to probe measurements filtered out. In other words, this is a version
%   of 'matches' without elements containing `NaN` values for any fields.
%
% matches -- All detected edges between coloured bands
%   A 2D structure array describing the matches found from detected edges
%   between bands of colour on the probe in the image to user-provided
%   measurements of the edges between probe colour bands. Each element
%   contains the match information for one detected edge. If no suitable
%   match was found for a given edge detected in the image, the
%   corresponding field values will be `NaN`. The fields of the structure
%   array are as follows:
%   - index: The index of the detected edge between two bands of colour
%       on the probe. The indices of edges correspond to the order of the
%       edges along the major axis of the probe in the image, from left to
%       right.
%   - lengthAlongPCAMajorAxis: The pixel distance along the major axis of
%       the probe from the first detected edge to the current edge. For
%       more information, refer to the documentation of the `lengths`
%       output argument of 'bilateralModel()'.
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
%       a structure input and output by 'CreateProbeDetectionModel.m'.
%   - matchedLength: The value of `probe.lengths(matchedLengthIndex)`, where
%       `probe` is a structure input and output by
%       'CreateProbeDetectionModel.m'.
%   - matchedWidth: The value of `probe.widths(matchedLengthIndex)`, where
%       `probe` is a structure input and output by
%       'CreateProbeDetectionModel.m'.
%
%   Each column of `matches` is a separate data association hypothesis - A
%   possible alignment of the detected edges to the edges of the probe.
%   Consequently, some columns may have elements in common. Data
%   association hypotheses are selected for having the largest number of
%   consistent pairings between detected and known edges, as computed by
%   'matchProbeLengthsRandom()'.
%
% ## Algorithm
%
% Probe detection begins with two colour detection passes: First, coloured
% regions which may correspond to probe bands are detected in the image,
% and bounding boxes are created to surround them. Alternatively, the user
% can provide a bounding region, as determined by
% `params.request_bounding_region`. Second, the same colours are detected
% within the bounding boxes, normally using less conservative parameters in
% order to avoid eroding the detected probe regions.
%
% Next, coloured regions are interpreted as bands along the length of the
% probe, and the endpoints of edges between bands are extracted from the
% borders of the coloured regions. The resulting sequence of measurements
% in the image is then matched to the user-provided measurements of the
% probe. At the end of this data association step, the probe could be
% overlayed on the image, but its location in 3D space is undetermined.
%
% ## Notes
% - If `params.request_bounding_region` is `true`, this function will open
%   a figure to collect a region of interest from the user. Therefore,
%   `params.request_bounding_region` should be `false` when this function
%   is intended to be non-interactive.
%
% See also roipoly, detectBoundingBoxes, detectWithinBoundingBox, bilateralModel, labelDetectedColors, matchProbeLengthsRandom

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created January 9, 2018

%% Parse input arguments

nargoutchk(1,2);
narginchk(5,6);

detectBoundingBoxesVerbose = [];
detectWithinBoundingBoxVerbose = [];
if ~isempty(varargin)
    verbose = varargin{1};
    if isfield(verbose, 'detectBoundingBoxesVerbose')
        detectBoundingBoxesVerbose = verbose.detectBoundingBoxesVerbose;
    end
    if isfield(verbose, 'detectWithinBoundingBoxVerbose')
        detectWithinBoundingBoxVerbose = verbose.detectWithinBoundingBoxVerbose;
    end
    display_detected_model_from_image = verbose.display_detected_model_from_image;
    display_final_clipped_regions_colored = verbose.display_final_clipped_regions_colored;
    verbose_detected_point_sequence_matching = verbose.verbose_detected_point_sequence_matching;
    display_detected_model_matching = verbose.display_detected_model_matching;
    warnings = verbose.warnings;
else
    display_detected_model_from_image = false;
    display_final_clipped_regions_colored = false;
    verbose_detected_point_sequence_matching = false;
    display_detected_model_matching = false;
    warnings = false;
end

if isempty(detectBoundingBoxesVerbose)
    detectBoundingBoxesVerbose.display_original_image = false;
    detectBoundingBoxesVerbose.display_hue_image = false;
    detectBoundingBoxesVerbose.display_saturation_image = false;
    detectBoundingBoxesVerbose.plot_hue_estimator = false;
    detectBoundingBoxesVerbose.plot_hue_classifier = false;
    detectBoundingBoxesVerbose.display_distribution_backprojections = false;
    detectBoundingBoxesVerbose.display_binary_images = false;
    detectBoundingBoxesVerbose.verbose_region_filtering = false;
    detectBoundingBoxesVerbose.display_region_expansion = false;
end

if isempty(detectWithinBoundingBoxVerbose)
    detectWithinBoundingBoxVerbose.display_hue_image = false;
    detectWithinBoundingBoxVerbose.display_saturation_image = false;
    detectWithinBoundingBoxVerbose.plot_hue_estimator = false;
    detectWithinBoundingBoxVerbose.plot_hue_classifier = false;
    detectWithinBoundingBoxVerbose.display_distribution_backprojections = false;
    detectWithinBoundingBoxVerbose.display_binary_images = false;
    detectWithinBoundingBoxVerbose.verbose_region_filtering = false;
    detectWithinBoundingBoxVerbose.display_regions_colored = false;
    detectWithinBoundingBoxVerbose.display_band_edge_extraction = false;
    detectWithinBoundingBoxVerbose.verbose_edge_endpoint_extraction = false;
end

image_width = size(I, 2);
image_height = size(I, 1);
image_n_channels = size(I, 3);
if image_n_channels ~= 3
    error('A 3-channel RGB image is required.')
end

%% Manual bounding region selection (if enabled)
if params.request_bounding_region
    figure
    imshow(I);
    title('Select a bounding region')
    probe_mask_initial = roipoly;
    if isempty(probe_mask_initial)
        error('Bounding region selection aborted. Either try again, or set `params.request_bounding_region` to false.');
    end
else
    if params.uniform_background_initial
        rgb_sigma_polyfit_initial = [];
    else
        rgb_sigma_polyfit_initial = rgb_sigma_polyfit;
    end
    probe_mask_initial = detectBoundingBoxes(...
        I, probe,...
        probe_color_distributions, rgb_sigma_polyfit_initial,...
        params.detectBoundingBoxesParams, detectBoundingBoxesVerbose...
    );
end

%% Detect interest points and refined probe colour regions within the bounding box

if params.uniform_background_final
    rgb_sigma_polyfit_final = [];
else
    rgb_sigma_polyfit_final = rgb_sigma_polyfit;
end
[...
    interest_points_detected, probe_regions_bw_final_filtered...
] = detectWithinBoundingBox(...
    I, probe, probe_mask_initial,...
    probe_color_distributions, rgb_sigma_polyfit_final,...
    params.detectWithinBoundingBoxParams, detectWithinBoundingBoxVerbose...
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
    interest_points_detected, params.detected_point_alignment_outlier_threshold, false...
);

n_detected_band_edges = length(image_lengths_detected);

if display_detected_model_from_image
    fg = figure;
    imshow(I);
    plotBilateralModel( model_from_image_detected, model_from_image_detected_axes, [image_height, image_width], [], fg );
    title('Classified detected interest points (green, red = above/below first PCA component; yellow = unmatched)');
end

if warnings && ~isempty(model_from_image_detected.unmatched)
    warning(sprintf(['Not all detected probe colour band junctions are marked with two points - One on each edge of the probe.\n',...
        'The output of probe colour band junction detection may have been quite messy.\n',...
        'Consider adjusting edge detection and edge endpoint detection parameters,\n',...
        'or increasing the outlier detection threshold used when pairing points.'])); %#ok<SPWRN>
end

%% Label the colours between detected probe band edges

[...
    probe_colors_detected_left, probe_colors_detected_right...
] = labelDetectedColors(...
    I, probe_color_distributions,...
    probe_regions_bw_final_filtered,...
    model_from_image_detected_pca_space, model_from_image_detected,...
    model_from_image_transform, model_from_image_detected_axes(1, :),...
    display_final_clipped_regions_colored...
);

%% Match model extracted from the image to the user-supplied measurements of the probe

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

alignment_thresholds.direction_threshold = params.matching.direction_threshold;
alignment_thresholds.inlier_threshold = params.matching.inlier_threshold;
alignment_thresholds.local_alignment_threshold = params.matching.local_alignment_threshold;
image_to_measured_matches_detected = matchProbeLengthsRandom(...
  subject(:, 1),...
  query(:, 1),...
  params.matching.subject_gap_cost_detection,...
  params.matching.query_gap_cost_detection,...
  subject(:, 2:3),...
  query(:, 2:3),...
  alignment_thresholds,...
  verbose_detected_point_sequence_matching...
);

% Remove probe tips
tip_filter = image_to_measured_matches_detected == matching_start_index;
if any(any(tip_filter))
    if warnings
        warning('Removing detected interest points which were matched to the start of the probe.')
    end
    image_to_measured_matches_detected(tip_filter) = 0;
end
tip_filter = image_to_measured_matches_detected == matching_end_index;
if any(any(tip_filter))
    if warnings
        warning('Removing detected interest points which were matched to the end of the probe.')
    end
    image_to_measured_matches_detected(tip_filter) = 0;
end

%% Output

probe_widths_for_matching = probe.widths(matching_start_index:matching_end_index);
image_to_measured_matches_detected_filter = logical(image_to_measured_matches_detected);
image_to_measured_matches_detected(~image_to_measured_matches_detected_filter) = NaN;
n_hypotheses = size(image_to_measured_matches_detected_filter, 2);
if ~all(all(image_to_measured_matches_detected_filter))
    if warnings
        warning('Some detected interest points in the image were not matched to known probe measurements.')
    end
    matched_lengths = nan(n_detected_band_edges, n_hypotheses);
    matched_lengths(image_to_measured_matches_detected_filter) =...
        probe_lengths_for_matching(...
                image_to_measured_matches_detected(...
                        image_to_measured_matches_detected_filter...
                    )...
            );
    matched_widths = nan(n_detected_band_edges, n_hypotheses);
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

matches = struct(...
        'index', repmat(num2cell((1:n_detected_band_edges).'), 1, n_hypotheses),...
        'lengthAlongPCAMajorAxis', repmat(num2cell(image_lengths_detected), 1, n_hypotheses),...
        'pointAbovePCAMajorAxis', repmat(num2cell(model_from_image_detected.above, 2), 1, n_hypotheses),...
        'pointBelowPCAMajorAxis', repmat(num2cell(model_from_image_detected.below, 2), 1, n_hypotheses),...
        'matchedLengthIndex', num2cell(...
                image_to_measured_matches_detected + (matching_start_index - 1)...
            ),...
        'matchedLength', num2cell(matched_lengths),...
        'matchedWidth', num2cell(matched_widths)...
    );

if display_detected_model_matching
    disp('Match between detected interest points and probe measurements:')
    probe_detection_matches_display =...
        rmfield(matches, 'lengthAlongPCAMajorAxis');
    probe_detection_matches_display =...
        rmfield(probe_detection_matches_display, 'matchedWidth');
    for hyp = 1:n_hypotheses
        fprintf('Hypothesis: %d\n', hyp)
        probe_detection_matches_display_hyp = struct2table(probe_detection_matches_display(:, hyp));
        disp(probe_detection_matches_display_hyp);
    end
end

matches_filtered = cell(n_hypotheses, 1);
for hyp = 1:n_hypotheses
    matches_filtered{hyp} = matches(image_to_measured_matches_detected_filter(:, hyp), hyp);
end