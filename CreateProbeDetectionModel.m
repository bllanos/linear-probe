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

% Bernard Llanos
% Spring 2016 research assistantship supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created July 26, 2016

%% Input data and parameters

% Probe measurements
model_filename = 'C:\Users\Bernard\Documents\Data\20160725_probes_lotus\redPen.mat';
% Image of probe
I_filename = 'C:\Users\Bernard\Documents\Data\20160725_probes_lotus\probesAgainstWhiteBox\original\redPen_white_b_1.bmp';
% Annotations for image of probe
I_annotations_filename = 'C:\Users\Bernard\Documents\Data\20160725_probes_lotus\probesAgainstWhiteBox\annotated\redPen_white_b_1.png';

% Annotation extraction parameters
annotation_corner_search_width = 4; % Set to zero to use centers of user-marked annotations as opposed to nearby corner features

% Parameters for interpreting annotated points
point_alignment_outlier_threshold = 3;

% Parameters for matching annotated points with the probe measurements
subject_gap_cost = -0.1;
query_gap_cost = 0;

% Debugging tools
display_original_image = false;
display_annotations_image = false;
display_extracted_annotations = false;
display_model_from_image = false;
verbose_point_sequence_matching = false;

%% Load images and obtain adjusted centers of user-marked annotations

I = imread(I_filename);
image_width = size(I, 2);
image_height = size(I, 1);
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

% The probe tips are marked with single points. All other probe segments
% must be marked with a point along each side of the probe.
%
% Since we are performing calibration, failure to satisfy this condition is
% an error. In a detection scenario, we would just eliminate the outlier
% points and continue.
if ~isempty(model_from_image.unmatched)
    error('Not all probe color band junctions are marked with two points - One on each edge of the probe.');
end

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

%% Match model extracted from the image to the user-supplied measurements of the probe
load(model_filename);
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
  verbose_point_sequence_matching...
);