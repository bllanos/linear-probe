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
%     length of the entire probe).
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

% Debugging tools
display_original_image = false;
display_annotations_image = false;
display_extracted_annotations = false;

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
        title('Extracted user-marked interest points (red o) and other corner features');
        for rect = feature_search_windows'
            rectangle('Position', rect, 'EdgeColor', 'b');
        end
        for i = 1:length(corners)
            if ~isempty(corners{i})
                plot(corners{i});
            end
        end
    else
        title('Extracted user-marked interest points (red o)');
    end
    hold off
end




%% Model the interest points as a series of segments

% Determine probe centerline with linear regression
X = [ones(n_annotations, 1), feature_points(:, 1)];
intercept_slope = X \ feature_points(:, 2);
hold on
y_line = X * intercept_slope;
plot(feature_points(:, 1),y_line,'g--')
hold off

% Express interest points in a coordinate space where the centerline is the
% x-axis, and the leftmost interest point's closest projection onto the
% line is on the y-axis
% See https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line
a = intercept_slope(2);
b = -1;
c = intercept_slope(1);
feature_points_on_line_x = (b * (b * feature_points(:, 1) - a * feature_points(:, 2)) - a * c)...
    / (a^2 + b^2);
feature_points_on_line_y = (a * (-b * feature_points(:, 1) + a * feature_points(:, 2)) - b * c)...
    / (a^2 + b^2);
% distances_to_line = sqrt((feature_points(:, 1) - feature_points_on_line_x) .^ 2 ...
%     + (feature_points(:, 2) - feature_points_on_line_y) .^ 2);
[x_min_line, min_ind] = min(feature_points_on_line_x);
line_origin = [x_min_line, feature_points_on_line_y(min_ind)];
line_direction = [1, intercept_slope(2)];
line_direction = line_direction / norm(line_direction);
line_direction_rep = repmat(line_direction, n_annotations, 1);
feature_points_to_line_origin = feature_points - repmat(line_origin, n_annotations, 1);
feature_points_line_coordinates_x = dot(feature_points_to_line_origin, line_direction_rep, 2);
feature_points_line_coordinates_y = feature_points_to_line_origin - repmat(feature_points_line_coordinates_x, 1, 2) .* line_direction_rep;
feature_points_line_coordinates_y = sign(feature_points_line_coordinates_y(:, 2))...
    .* sqrt(dot(feature_points_line_coordinates_y, feature_points_line_coordinates_y, 2));

% Partition points into those above and below the line
above_line_filter = feature_points_line_coordinates_y > 0;
feature_points_above = feature_points(above_line_filter, :);
feature_points_below = feature_points(~above_line_filter, :);
feature_points_line_coordinates_above = [
        feature_points_line_coordinates_x(above_line_filter),...
        feature_points_line_coordinates_y(above_line_filter)
    ];
feature_points_line_coordinates_below = [
        feature_points_line_coordinates_x(~above_line_filter),...
        feature_points_line_coordinates_y(~above_line_filter)
    ];

% Sort points by x-coordinate
A = [feature_points_line_coordinates_above, feature_points_above];
A = sortrows(A);
feature_points_line_coordinates_above = A(:, 1:2);
feature_points_above = A(:, 3:4);

A = [feature_points_line_coordinates_below, feature_points_below];
A = sortrows(A);
feature_points_line_coordinates_below = A(:, 1:2);
feature_points_below = A(:, 3:4);

figure
hold on
plot(feature_points_line_coordinates_x, zeros(n_annotations, 1),'k--')
plot(feature_points_line_coordinates_above(:, 1), feature_points_line_coordinates_above(:, 2), 'g.')
plot(feature_points_line_coordinates_below(:, 1), feature_points_line_coordinates_below(:, 2), 'r.')
title('Probe interest points in centerline coordinate space')
xlabel('X')
ylabel('Y')

% Collapse the points into a single row by matching points above and below
% the line
n_feature_points_above = size(feature_points_line_coordinates_above, 1);
min_x_separations_above = zeros(n_feature_points_above, 1);
min_x_separation_indices_above = zeros(n_feature_points_above, 1);
for i = 1:n_feature_points_above
    [min_x_separations_above(i), min_x_separation_indices_above(i)] =...
        min(abs(feature_points_line_coordinates_below(:, 1) - feature_points_line_coordinates_above(i, 1)));
end

n_feature_points_below = size(feature_points_line_coordinates_below, 1);
min_x_separations_below = zeros(n_feature_points_below, 1);
min_x_separation_indices_below = zeros(n_feature_points_below, 1);
for i = 1:n_feature_points_below
    [min_x_separations_below(i), min_x_separation_indices_below(i)] =...
        min(abs(feature_points_line_coordinates_above(:, 1) - feature_points_line_coordinates_below(i, 1)));
end

% Identify mismatched points based on mutual selection
inliers_above = (min_x_separation_indices_below(min_x_separation_indices_above) == (1:n_feature_points_above).');
inliers_below = (min_x_separation_indices_above(min_x_separation_indices_below) == (1:n_feature_points_below).');

% Identify mismatched points with a statistical analysis
% Reference: MATLAB documentation page on "Inconsistent Data"

% Outlier separation in x-coordinates
inlier_separations = [min_x_separations_above(inliers_above); min_x_separations_below(inliers_below)];
mu = mean(inlier_separations);
sigma = std(inlier_separations);

mu_rep = repmat(mu,n_feature_points_above,1);
sigma_rep = repmat(sigma,n_feature_points_above,1);
point_alignment_outlier_threshold = 3;
fprintf('(Parameter) Outlier threshold (number of standard deviations) for detecting interest points on different bands: %g\n', point_alignment_outlier_threshold);
inliers_above = inliers_above & (abs(min_x_separations_above - mu_rep) < point_alignment_outlier_threshold * sigma_rep);

mu_rep = repmat(mu,n_feature_points_below,1);
sigma_rep = repmat(sigma,n_feature_points_below,1);
inliers_below = inliers_below & (abs(min_x_separations_below - mu_rep) < point_alignment_outlier_threshold * sigma_rep);

% Outlier separation in y-coordinates
y_separations_above = abs(feature_points_line_coordinates_above(:, 2)...
    + feature_points_line_coordinates_below(min_x_separation_indices_above, 2));
y_separations_below = abs(feature_points_line_coordinates_below(:, 2)...
    + feature_points_line_coordinates_above(min_x_separation_indices_below, 2));
inlier_separations = abs([
    feature_points_line_coordinates_above(inliers_above, 2)...
    + feature_points_line_coordinates_below(min_x_separation_indices_above(inliers_above), 2);
    feature_points_line_coordinates_below(inliers_below, 2)...
    + feature_points_line_coordinates_above(min_x_separation_indices_below(inliers_below), 2);]);
mu = mean(inlier_separations);
sigma = std(inlier_separations);

mu_rep = repmat(mu,n_feature_points_above,1);
sigma_rep = repmat(sigma,n_feature_points_above,1);
inliers_above = inliers_above & (abs(y_separations_above - mu_rep) < point_alignment_outlier_threshold * sigma_rep);

mu_rep = repmat(mu,n_feature_points_below,1);
sigma_rep = repmat(sigma,n_feature_points_below,1);
inliers_below = inliers_below & (abs(y_separations_below - mu_rep) < point_alignment_outlier_threshold * sigma_rep);

% It is acceptable if the probe tips are marked with single points.
% All other probe segments must be marked with a point on each side of the
% centerline.
%
% Since we are performing calibration, failure to satisfy this condition is
% an error. In a detection scenario, we would just eliminate the outlier
% points and continue.
A = [
    feature_points_line_coordinates_above(:, 1), double(inliers_above);
    feature_points_line_coordinates_below(:, 1), double(inliers_below)
    ];
A = sortrows(A);
all_inliers = A(:, end);
n_outliers = sum(~all_inliers);
if (~all_inliers(1) + ~all_inliers(end)) ~= n_outliers
    error('Not all probe segment junctions are marked with two points, above and below the centerline.');
end
probe_has_left_tip = ~all_inliers(1);
probe_has_right_tip = ~all_inliers(end);

% Collect the single row of points
n_boundary_points = n_feature_points_above + sum(~inliers_below);
boundary_points = zeros(n_boundary_points, 2);
segments_to_points_above = zeros(n_boundary_points, 1);
segments_to_points_below = zeros(n_boundary_points, 1);
start_index = 1;
index_above = 1;
end_index = n_boundary_points;
if probe_has_left_tip
    boundary_points(1, :) = 0;
    if ~inliers_above(1)
        segments_to_points_above(1) = 1;
        index_above = 2;
    else
        segments_to_points_below(1) = 1;
    end
    start_index = 2;
end
if probe_has_right_tip
    boundary_points(end, 2) = 0;
    if ~inliers_above(end)
        segments_to_points_above(end) = n_feature_points_above;
        boundary_points(end, 1) = feature_points_line_coordinates_above(end, 1);
    else
        segments_to_points_below(end) = n_feature_points_below;
        boundary_points(end, 1) = feature_points_line_coordinates_below(end, 1);
    end
    end_index = n_boundary_points - 1;
end

for i = start_index:end_index
    match_index = min_x_separation_indices_above(index_above);
    boundary_points(i, 1) = mean([
        feature_points_line_coordinates_above(index_above, 1), feature_points_line_coordinates_below(match_index, 1)
        ]);
    boundary_points(i, 2) = feature_points_line_coordinates_above(index_above, 2) - feature_points_line_coordinates_below(match_index, 2);
    index_above = index_above + 1;
end

% Plot for validation
plot(boundary_points(:, 1), boundary_points(:, 2), 'b.')
legend('Regression centerline', 'Interest points above line',...
    'Interest points below line', 'Merged interest points with double width');
hold off

%% Match segment model to the user-supplied model of the probe
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