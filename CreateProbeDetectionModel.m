%% Detecting a red pen probe - Rough work

% Bernard Llanos
% Spring 2016 research assistantship supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created July 26, 2016

clear
close all

I = imread('C:\Users\Bernard\Documents\Data\20160725_probes_lotus\probesAgainstWhiteBox\original\redPen_white_b_1.bmp');
image_width = size(I, 2);
image_height = size(I, 1);
[~, ~, Ialpha] = imread('C:\Users\Bernard\Documents\Data\20160725_probes_lotus\probesAgainstWhiteBox\annotated\redPen_white_b_1.png');
figure; imshow(I); title('Raw image');

%% Obtain centers of user-marked annotations

figure; imshow(Ialpha); title('User-marked interest points');
annotation = logical(Ialpha);
[annotations_y, annotations_x] = find(bwmorph(annotation, 'shrink', Inf));
hold on
plot(annotations_x, annotations_y, 'ro')
hold off

%% Test various corner detectors to try and refine annotations

n_annotations = length(annotations_x);
delta = 4;
width = delta * 2 + 1;
fprintf('(Parameter) Maximum feature point search window width: %d\n', width);

I_grey = rgb2gray(I);

annotation_search_windows = zeros(n_annotations, 4);
for i = 1:n_annotations
    x_separation = annotations_x([1:(i-1), (i+1):end]) - annotations_x(i);
    y_separation = annotations_y([1:(i-1), (i+1):end]) - annotations_y(i);
    distances = sqrt(x_separation .^ 2 + y_separation .^ 2);
    delta_i = min(floor((min(distances) - 1) / 2), delta);
    delta_i = max(0, delta_i);
    width_i = delta_i * 2 + 1;
    y_min = max(1, annotations_y(i) - delta_i);
    x_min = max(1, annotations_x(i) - delta_i);
    annotation_search_windows(i, :) = [x_min, y_min, width_i, width_i];
end

fast_corners = cell(n_annotations, 1);
min_eigen_corners = cell(n_annotations, 1);
harris_corners = cell(n_annotations, 1);
brisk_corners = cell(n_annotations, 1);
for i = 1:n_annotations
    fast_corners{i} = detectFASTFeatures(I_grey, 'ROI', annotation_search_windows(i, :));
    min_eigen_corners{i} = detectMinEigenFeatures(I_grey, 'ROI', annotation_search_windows(i, :));
    harris_corners{i} = detectHarrisFeatures(I_grey, 'ROI', annotation_search_windows(i, :));
    brisk_corners{i} = detectBRISKFeatures(I_grey, 'ROI', annotation_search_windows(i, :));
end

% figure; imshow(max(I_grey, Ialpha)); title('Feature detection results: FAST');
% hold on
% for i = 1:n_annotations
%     rectangle('Position', annotation_search_windows(i, :), 'EdgeColor', 'w');
%     if ~isempty(fast_corners{i})
%         plot(fast_corners{i});
%     end
% end
% hold off
% 
% figure; imshow(max(I_grey, Ialpha)); title('Feature detection results: minEigen');
% hold on
% for i = 1:n_annotations
%     rectangle('Position', annotation_search_windows(i, :), 'EdgeColor', 'w');
%     if ~isempty(min_eigen_corners{i})
%         plot(min_eigen_corners{i});
%     end
% end
% hold off
% 
% figure; imshow(max(I_grey, Ialpha)); title('Feature detection results: Harris');
% hold on
% for i = 1:n_annotations
%     rectangle('Position', annotation_search_windows(i, :), 'EdgeColor', 'w');
%     if ~isempty(harris_corners{i})
%         plot(harris_corners{i});
%     end
% end
% hold off
% 
% figure; imshow(max(I_grey, Ialpha)); title('Feature detection results: BRISK');
% hold on
% for i = 1:n_annotations
%     rectangle('Position', annotation_search_windows(i, :), 'EdgeColor', 'w');
%     if ~isempty(brisk_corners{i})
%         plot(brisk_corners{i});
%     end
% end
% hold off

% minEigen features seem good
% - Enough are detected
% - They are close to the user-specified points
feature_points = zeros(n_annotations, 2);
chosen_features = min_eigen_corners;
for i = 1:n_annotations
    if ~isempty(chosen_features{i})
        feature_points(i, :) = chosen_features{i}.selectStrongest(1).Location;
    else
        feature_points(i, :) = [annotations_x(i), annotations_y(i)];
    end
end

% Plot final point selections
figure; imshow(max(I_grey, Ialpha)); title('Adjusted interest points');
hold on
plot(feature_points(:, 1), feature_points(:, 2), 'g.')
for i = 1:n_annotations
    rectangle('Position', annotation_search_windows(i, :), 'EdgeColor', 'w');
end
hold off

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