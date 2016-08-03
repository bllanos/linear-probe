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
%   Specifically asymmetry of band lengths is assumed, as opposed to colour
%   pattern asymmetry. Band edges should be perpendicular to the probe's
%   axis of cylindrical symmetry.
%
% ## Input
%
% ### Probe measurements
% A '.mat' file containing a structure called 'probe' with the following
% fields:
% - lengths: Distances of edges of bands from the active end of the probe,
%     including a distance of 0.0 for the active end of the probe, and,
%     optionally, a distance for the other end of the probe (i.e. the
%     length of the entire probe). Distances are measured along the probe's
%     axis of cylindrical symmetry.
% - widths: Width of the probe at the edges of bands. Widths must include
%     the ends of the probe, except for ends that taper to points. The
%     elements of 'widths' should correspond to elements of 'lengths'
%     (although, if the probe tips taper to points, the first and/or last
%     elements of 'lengths' will not match any elements of 'widths').
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
% the probe's contour in the image.

% Bernard Llanos
% Spring 2016 research assistantship supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created July 26, 2016

%% Input data and parameters

% Image of probe
I_filename = 'C:\Users\Bernard\Documents\Data\20160725_probes_lotus\probesAgainstWhiteBox\original\redPen_white_b_1.bmp';
% Annotations for image of probe
I_annotations_filename = 'C:\Users\Bernard\Documents\Data\20160725_probes_lotus\probesAgainstWhiteBox\annotated\redPen_white_b_1.png';

% Annotation extraction parameters
annotation_corner_search_width = 4; % Set to zero to use centers of user-marked annotations as opposed to nearby corner features

% Parameters for interpreting annotated points
point_alignment_outlier_threshold = 3;

% Debugging tools
display_original_image = false;
display_annotations_image = false;
display_extracted_annotations = false;
display_model_from_image = true;

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
    figure;
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
load('C:\Users\Bernard\Documents\Data\20160725_probes_lotus\redPen.mat');
if ~exist('probe', 'var')
    error('No variable called ''probe'' is loaded.')
end

% The cross ratio is invariant under projective transformations
% See https://en.wikipedia.org/wiki/Cross-ratio
%
% Without knowing the pose of the probe in 3D space, it is only possible to
% compare the cross ratios of points along the probe with the cross ratios
% of points along the image of the probe. Lengths and length ratios cannot
% be compared directly.

% Computing the cross ratio requires 4 points
if length(probe.lengths) < 4
    error('Insufficient points given in the probe model to compute cross ratios.')
end
if n_boundary_points < 4
    error('Insufficient segments extracted from the image to compute cross ratios.')
end

% Compute all possible cross ratios of the image segments and probe points
% Note that 'nchoosek' preserves the order of the items being chosen.
image_combinations = nchoosek(1:n_boundary_points, 4);
n_image_cross_ratios = size(image_combinations, 1);
image_cross_ratios = zeros(n_image_cross_ratios, 1);
for i = 1:n_image_cross_ratios
    points = boundary_points(image_combinations(i, :), 1);
    image_cross_ratios(i) = crossRatio(points);
end

probe_lengths = probe.lengths;
n_probe_lengths = length(probe_lengths);
probe_combinations = nchoosek(1:n_probe_lengths, 4);
n_probe_cross_ratios = size(probe_combinations, 1);
probe_cross_ratios = zeros(n_probe_cross_ratios, 1);
for i = 1:n_probe_cross_ratios
    points = probe_lengths(probe_combinations(i, :), 1);
    probe_cross_ratios(i) = crossRatio(points);
end

% Plot cross ratios
figure;
hold on
plot(image_cross_ratios, 'r.');
plot(probe_cross_ratios, 'g.');
hold off
title('Cross ratios')
xlabel('Combination index')
ylabel('Cross ratio')
legend('Image cross ratios', 'Measured probe cross ratios')

% Group cross ratios into closest pairs
min_cross_ratio_image_to_probe = zeros(n_image_cross_ratios, 1);
min_cross_ratio_image_to_probe_indices = zeros(n_image_cross_ratios, 1);
for i = 1:n_image_cross_ratios
    [min_cross_ratio_image_to_probe(i), min_cross_ratio_image_to_probe_indices(i)] =...
        min(abs(probe_cross_ratios - image_cross_ratios(i)));
end

min_cross_ratio_probe_to_image = zeros(n_probe_cross_ratios, 1);
min_cross_ratio_probe_to_image_indices = zeros(n_probe_cross_ratios, 1);
for i = 1:n_probe_cross_ratios
    [min_cross_ratio_probe_to_image(i), min_cross_ratio_probe_to_image_indices(i)] =...
        min(abs(image_cross_ratios - probe_cross_ratios(i)));
end

% Identify potentially mismatched cross ratios based on mutual selection
% This will filter out a lot of valid matches due to the large number of 
inliers_image_to_probe_filter = (...
        min_cross_ratio_probe_to_image_indices(min_cross_ratio_image_to_probe_indices)...
        == (1:n_image_cross_ratios).' ...
    );

% Final matching of cross ratios
inlier_image_combinations = image_combinations(inliers_image_to_probe_filter, :);
inlier_probe_combinations = probe_combinations(min_cross_ratio_image_to_probe_indices(inliers_image_to_probe_filter), :);

% For each pairing of cross ratios, increment the vote for a match between
% the corresponding probe and image points
%
% This needs to be done twice, because the cross ratio is the same for
% points in reverse order.
image_to_probe_match_votes = zeros(n_boundary_points, n_probe_lengths);
reversed_image_to_probe_match_votes = zeros(n_boundary_points, n_probe_lengths);
reversed_inlier_image_combinations = fliplr(inlier_image_combinations);