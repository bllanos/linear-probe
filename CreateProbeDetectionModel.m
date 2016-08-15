%% Probe Object Detection Model Creation
% Associate user-provided measurements of the probe with points marked on
% the image of the probe, such that the image colours can be used to detect
% the probe in other images.
%
% ## Usage
%   Modify parameters and paths to input data in the first code section
%   below, then run.
%
% ## Probe Design Assumptions
%   The probe is a cylindrically-symmetric object composed of conical or
%   cylindrical segments. It may taper to a point at one or both ends,
%   although only one end will be used as the active tip.
%
%   The length of the probe should be divided into different-coloured
%   bands, and the pattern of bands should be asymmetrical, such that it is
%   possible to uniquely determine the orientation of the probe.
%   Specifically, asymmetry of band lengths is assumed, as opposed to colour
%   pattern asymmetry. Band edges should be perpendicular to the probe's
%   axis of cylindrical symmetry.
%
% ## Input
%
% ### Probe measurements
% A '.mat' file containing a structure called 'probe' with the following
% fields:
% - lengths: Distances of edges of bands from the active end of the probe,
%     including a distance of 0.0 for the active end of the probe, and a
%     distance for the other end of the probe (i.e. the length of the
%     entire probe). Distances are measured along the probe's axis of
%     cylindrical symmetry.
% - colors: A vector with `length(lengths) - 1` elements specifying colour
%   indices for the bands of the probe. Colour indices should be
%   consecutive integers starting at 1. They indicate how the bands of
%   the probe are grouped based on mutually-distinguishable colours,
%   allowing bands to have non-unique colours. The specific index assigned
%   to a given band is unimportant.
% - widths: Width of the probe at the edges of bands. Widths must include
%     the ends of the probe, with values of zero for ends that taper to
%     points. The elements of 'widths' should correspond to elements of
%     'lengths' as follows:
%     1) Active tip of probe.
%     2) `2 * (length(probe.lengths) - 2)` measurements: One for each side
%        of each edge between bands, to disambiguate between conical and
%        cylindrical sections and to allow for junctions with slightly different
%        radii on either side. Measurements are in order of distance from
%        the active tip of the probe.
%     3) Other end of the probe.
%
% Units are arbitrary, but should be consistent with the units of any
% partial reconstruction of an object that the probe is used to refine.
%
% ### Image of probe
% An RGB image in any format that can be loaded by `imread` showing the
% entire portion of the probe that the user has provided measurements for.
% The image should depict the probe under lighting conditions that are
% similar to the detection scenario, to facilitate colour-based detection.
%
% ### Annotations for probe image
% An image in any format that can be loaded by `imread` with an alpha
% channel that is nonzero at user-marked interest points, and zero
% everywhere else. The image should be based on the same image of the probe
% provided above.
%
% Interest points can be marked by nonzero alpha channel regions that are
% more than one pixel in size - A morphological shrink operation will be
% used to extract single pixel locations from them. Optionally, a corner
% feature detector will then be used to refine the locations of interest
% points within a search window.
%
% A single interest point should be marked for each end of the probe that
% tapers to a point. Two interest points should be marked for the edges of
% coloured bands on the probe, corresponding to their intersections with
% the probe's contour in the image. For junctions between coloured bands of
% different radii, two interest points should be marked on the same side
% of the junction. Assuming the difference in radii is not significant, the
% choice of side should not matter with respect to the determination of the
% average colour on each side of the junction.
%
% ### Colour noise parameters
% A '.mat' file containing a 'rgb_sigma_polyfit' variable, as output by the
% script '.\EstimateRGBStandardDeviations.m'. 'rgb_sigma_polyfit' describes
% the variation in RGB channel standard deviations with RGB values in the
% image. This information should be computed from images taken under the
% same conditions and with the same camera parameters as the image used for
% probe modeling, if not computed from the same image used for probe
% modeling.
%
% ## Output
%
% ### Probe detection model
% A '.mat' file containing the following variables:
% - 'probe': A copy of the 'probe' variable loaded from the probe
%   measurements file, for convenience.
% - 'probe_band_color_distributions': Discretized variable kernel density
%   estimators of image hue values corresponding to the different coloured
%   bands on the probe, in the same order (starting from the active tip of
%   the probe). The i-th column of this 2D array stores the estimator for
%   the i-th probe segment.
% - 'probe_band_color_distribution_increment': A scalar equal to the spacing
%   between the samples of hue values in the range [0, 1] at which the
%   variable kernel density estimators have been evaluated to produce
%   'probe_band_color_distributions'. To find the approximate value of the
%   i-th estimator at a query value 'x' in the range [0, 1], use:
%
%   ```
%     queryDiscretized1DFunction(...
%       x, probe_band_color_distributions(:, i),...
%       probe_band_color_distribution_increment...
%     )
%   ```
%
% Additionally, the output file contains the values of all parameters in
% the first section of the script below, for reference. (Specifically,
% those listed in `parameters_list`, which should be updated if the set of
% parameters is changed.)
%
% ## References
% - T. Gevers and H. Stokman. "Robust Histogram Construction from Color
%   Invariants for Object Recognition". IEEE Transactions on Pattern
%   Analysis and Machine Intelligence, vol. 26, no. 1, pp. 113-118, Jan.
%   2004.

% Bernard Llanos
% Spring 2016 research assistantship supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created July 26, 2016

%% Input data and parameters

% List of parameters to save with results
parameters_list = {
        'model_filename',...
        'I_filename',...
        'I_annotations_filename',...
        'rgb_sigma_filename',...
        'annotation_corner_search_width',...
        'point_alignment_outlier_threshold',...
        'subject_gap_cost',...
        'query_gap_cost',...
        'n_samples_sequence_alignment',...
        'probe_color_distribution_resolution'
    };

% Probe measurements
model_filename = 'C:\Users\Bernard\Documents\Data\20160811_bambooSkewerProbe\bambooSkewer_orangeBlue.mat';
% Image of probe
I_filename = 'C:\Users\Bernard\Documents\Data\20160811_bambooSkewerProbe\original\probePrePaperOcclusion_1_b.bmp';
% Annotations for image of probe
I_annotations_filename = 'C:\Users\Bernard\Documents\Data\20160811_bambooSkewerProbe\annotated\probePrePaperOcclusion_1_b.png';
% RGB noise parameters
rgb_sigma_filename = 'C:\Users\Bernard\Documents\Data\20160811_bambooSkewerProbe\20160811_rgbStddev_bottomCamera.mat';

% Annotation extraction parameters
annotation_corner_search_width = 4; % Set to zero to use centers of user-marked annotations as opposed to nearby corner features

% Parameters for interpreting annotated points
point_alignment_outlier_threshold = 5;

% Parameters for matching annotated points with the probe measurements
subject_gap_cost = -0.1;
query_gap_cost = 0;
n_samples_sequence_alignment = 8;

% Number of points at which to evaluate hue variable kernel density estimators
probe_color_distribution_resolution = 180;

% Debugging tools
display_original_image = false;
display_annotations_image = false;
display_extracted_annotations = false;
display_model_from_image = false;
verbose_point_sequence_matching = false;
display_probe_band_masks = false;
display_probe_color_masks = false;
display_hue_image = false;
plot_hue_estimators = false;

%% Load images and obtain adjusted centers of user-marked annotations

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

[~, ~, Ialpha] = imread(I_annotations_filename);
if display_annotations_image
    figure; %#ok<UNRCH>
    imshow(Ialpha);
    title('User-marked annotations');
end

if annotation_corner_search_width
    [ interest_points, feature_search_windows, corners ] = extractInterestPoints( Ialpha, I, annotation_corner_search_width );
else
    interest_points = extractInterestPoints( Ialpha ); %#ok<UNRCH>
end
if display_extracted_annotations
    figure; %#ok<UNRCH>
    I_grey = rgb2gray(I);
    imshow(max(I_grey, Ialpha));
    hold on
    plot(interest_points(:, 1), interest_points(:, 2), 'ro')
    if annotation_corner_search_width
        title('Final marked interest points (red o) and other corner features');
        for rect = feature_search_windows'
            rectangle('Position', rect, 'EdgeColor', 'b');
        end
        for i = 1:length(corners)
            if ~isempty(corners{i})
                plot(corners{i});
            end
        end
    else
        title('Extracted interest points (red o)');
    end
    hold off
end

%% Model the interest points as a series of bands on a linear probe object

[...
    model_from_image, model_from_image_lengths, model_from_image_axes, model_to_image_transform...
] = bilateralModel(...
    interest_points, point_alignment_outlier_threshold, true...
);

if display_model_from_image
    figure; %#ok<UNRCH>
    imshow(I);
    title('Classified interest points (blue, black = tips; green, red = above/below first PCA component; yellow = unmatched)');
    hold on
    
    % Plot points
    if isfield(model_from_image, 'head')
        head = [model_from_image.head, 1] * model_to_image_transform;
        scatter(head(1), head(2), 'b')
    end
    if isfield(model_from_image, 'tail')
        tail = [model_from_image.tail, 1] * model_to_image_transform;
        scatter(tail(1), tail(2), 'k')
    end
    above = [model_from_image.above, ones(size(model_from_image.above, 1), 1)] * model_to_image_transform;
    scatter(above(:, 1), above(:, 2), 'g')
    below = [model_from_image.below, ones(size(model_from_image.below, 1), 1)] * model_to_image_transform;
    scatter(below(:, 1), below(:, 2), 'r')
    unmatched = [model_from_image.unmatched, ones(size(model_from_image.unmatched, 1), 1)] * model_to_image_transform;
    scatter(unmatched(:, 1), unmatched(:, 2), 'y')
    
    % Plot PCA lines
    line_points = lineToBorderPoints(model_from_image_axes, [image_height, image_width]);
    line(line_points(1, [1,3])', line_points(1, [2,4])', 'Color', 'c');
    line(line_points(2, [1,3])', line_points(2, [2,4])', 'Color', 'm');
    hold off
end

% The probe tips are marked with single points. All other probe segments
% must be marked with a point along each side of the probe.
%
% Since we are performing calibration, failure to satisfy this condition is
% an error. In a detection scenario, we would just eliminate the outlier
% points and continue.
if ~isempty(model_from_image.unmatched)
    error(sprintf(['Not all probe color band junctions are marked with two points - One on each edge of the probe.\n',...
        'Consider increasing the outlier detection threshold used when pairing points.'])); %#ok<SPERR>
end

%% Match model extracted from the image to the user-supplied measurements of the probe
load(model_filename, 'probe');
if ~exist('probe', 'var')
    error('No variable called ''probe'' is loaded (which would contain probe measurements).')
end

image_lengths = mean([model_from_image.above(:, 1), model_from_image.below(:, 1)], 2);
if isfield(model_from_image, 'head')
    image_lengths = [model_from_image.head(1); image_lengths];
end
if isfield(model_from_image, 'tail')
    image_lengths = [image_lengths; model_from_image.tail(1)];
end
    
image_to_measured_matches = matchPointsByCrossRatios(...
  probe.lengths, image_lengths, subject_gap_cost, query_gap_cost,...
  n_samples_sequence_alignment, verbose_point_sequence_matching...
);

% Validate the matching, and flip the model extracted from the image if necessary
if any(~image_to_measured_matches)
    error('Some interest points in the image were not matched to known probe measurements.')
else
    image_to_measured_matches_expected = (1:length(image_to_measured_matches)).';
    aligned_forward = all(image_to_measured_matches == image_to_measured_matches_expected);
    aligned_reverse = all(image_to_measured_matches == flipud(image_to_measured_matches_expected));
    if ~xor(aligned_forward, aligned_reverse)
        error('Failed to find an ordered one-to-one mapping between interest points in the image and known probe measurements.')
    elseif aligned_reverse
        model_from_image_new = model_from_image;
        if isfield(model_from_image, 'head')
            model_from_image_new.tail = model_from_image.head;
        else
            model_from_image_new = rmfield(model_from_image_new, 'tail');
        end
        if isfield(model_from_image, 'tail')
            model_from_image_new.head = model_from_image.tail;
        else
            model_from_image_new = rmfield(model_from_image_new, 'head');
        end
        model_from_image_new.above = flipud(model_from_image_new.above);
        model_from_image_new.below = flipud(model_from_image_new.below);
        model_from_image = model_from_image_new;
    end
end

%% Locate regions corresponding to probe bands and probe colours

% Express the model in pixel coordinates as polygonal sections
n_bands = length(probe.lengths) - 1;
above = [model_from_image.above, ones(size(model_from_image.above, 1), 1)] * model_to_image_transform;
below = [model_from_image.below, ones(size(model_from_image.below, 1), 1)] * model_to_image_transform;
model_from_image_polygons = cell(n_bands, 1);
start_offset = 0;
if isfield(model_from_image, 'head')
    head = [model_from_image.head, 1] * model_to_image_transform;
    model_from_image_polygons{1} = [
            head(1:2);
            above(1, 1:2);
            below(1, 1:2)
        ];
    start_offset = 1;
end
for i = 1:(size(above, 1) - 1)
    model_from_image_polygons{i + start_offset} = [
            above(i:(i+1), 1:2);
            below((i+1):-1:i, 1:2)
        ];
end
if isfield(model_from_image, 'tail')
    tail = [model_from_image.tail, 1] * model_to_image_transform;
    model_from_image_polygons{end} = [
        above(end, 1:2);
        tail(1:2);
        below(end, 1:2);
    ];
end

% Obtain a mask for each band
probe_band_masks = false(image_height, image_width, n_bands);
for i = 1:n_bands
    probe_band_masks(:, :, i) = roipoly(...
            I, model_from_image_polygons{i}(:, 1), model_from_image_polygons{i}(:, 2)...
        );
end

if display_probe_band_masks
    for i = 1:n_bands %#ok<UNRCH>
        figure
        imshow(probe_band_masks(:, :, i));
        title(sprintf('Mask for probe band %d', i))
    end
end

% Group bands by colour
probe_colors = unique(probe.colors);
n_colors = length(probe_colors);
probe_color_masks = false(image_height, image_width, n_colors);
for i = 1:n_bands
    probe_color_masks(:, :, probe.colors(i)) =...
        probe_color_masks(:, :, probe.colors(i)) | probe_band_masks(:, :, i);
end

if display_probe_color_masks
    for i = 1:n_colors %#ok<UNRCH>
        figure
        imshow(probe_color_masks(:, :, i));
        title(sprintf('Mask for probe colour %d', i))
    end
end

%% Create photometric invariant representations of the probe colours

% Obtain hue values
H = rgb2hue(I);
    
if display_hue_image
    figure %#ok<UNRCH>
    imshow(H);
    title('Hue channel of original image')
end

% Compute hue variable kernel density estimators from probe colors
load(rgb_sigma_filename, 'rgb_sigma_polyfit');
if ~exist('rgb_sigma_polyfit', 'var')
    error('No variable called ''rgb_sigma_polyfit'' is loaded (which would contain the camera RGB noise model).')
end

probe_color_distributions = zeros(...
        probe_color_distribution_resolution, n_colors...
    );
I_double = im2double(I);
R = I_double(:, :, 1);
G = I_double(:, :, 2);
B = I_double(:, :, 3);
for i = 1:n_colors
    [...
        probe_color_distributions(:, i),...
        probe_color_distribution_increment...
    ] = hueVariableKernelDensityEstimator(...
        H, R, G, B, probe_color_masks(:, :, i),...
        rgb_sigma_polyfit, probe_color_distribution_resolution...
    );
end

if plot_hue_estimators
    x = 0:probe_color_distribution_increment:1; %#ok<UNRCH>
    line_styles = {'-', '--', ':', '-.'};
    legend_names = cell(n_colors, 1);
    plot_colors = jet(n_colors);
    figure
    hold on
    for i = 1:n_colors
        legend_names{i} = sprintf('Probe color %d', i);
        plot(...
                x, probe_color_distributions(:, i),...
                'Color', plot_colors(i, :),...
                'LineStyle', line_styles{mod(i - 1, length(line_styles)) + 1}...
            )
    end
    hold off
    legend(legend_names{:});
    title('Hue variable kernel density estimators for colors on the probe')
    xlabel('Hue, \theta (range [0, 1])')
    ylabel('Density, P(\theta)')
end

%% Save results to a file
save_variables_list = [ parameters_list, {...
        'probe',...
        'probe_color_distributions',...
        'probe_color_distribution_increment'...
    } ];
uisave(save_variables_list,'probeDetectionModel')
disp('Reminder: The output model is specific to the probe, camera, and camera parameters.')