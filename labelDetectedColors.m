function [ colors_left, colors_right ] = labelDetectedColors(...
    I, probe_color_distributions, probe_color_distribution_increment,...
    probe_regions_bw,...
    model_pca_space, model, model_transform, params, varargin...
)
%LABELDETECTEDCOLORS Robustly associate detected probe regions with colour classes
%
% ## Syntax
% [ colors_left, colors_right ] = labelDetectedColors(...
%     I, probe_color_distributions, probe_color_distribution_increment,...
%     probe_regions_bw,...
%     model_pca_space, model, model_transform, params [, verbose]...
% )
%
% ## Description
% [ colors_left, colors_right ] = labelDetectedColors(...
%     I, probe_color_distributions, probe_color_distribution_increment,...
%     probe_regions_bw,...
%     model_pca_space, model, model_transform, params [, verbose]...
% )
%   Returns the colour labels associated with the regions to either side of
%   each detected edge between the coloured bands of the probe.
%
% ## Input Arguments
%
% I -- Image containing the probe
%   An RGB image containing the probe.
%
% probe_color_distributions -- Probe colour estimators
%   Discretized variable kernel density estimators of image hue values
%   corresponding to the different coloured bands on the probe, in the same
%   order (starting from the active tip of the probe). The i-th column of
%   this 2D array stores the estimator for the i-th colour class of probe
%   segments. (There are `n_colors` columns in total.)
%
%   Presently, this argument is only used for producing graphical debugging
%   output, so an empty array could be passed if graphical debugging output
%   is disabled.
%
% probe_color_distribution_increment -- Probe colour estimator sample spacing
%   A scalar equal to the spacing between the samples of hue values in the
%   range [0, 1] at which the variable kernel density estimators have been
%   evaluated to produce 'probe_band_color_distributions'.
%
%   Presently, this argument is only used for producing graphical debugging
%   output, so an empty array could be passed if graphical debugging output
%   is disabled.
%
% probe_regions_bw -- Detected probe colour bands
%   A three-dimensional logical array of size image_height x image_width x
%   n_colors. `probe_regions_bw_filtered(:, :, i)` stores the binary image
%   representing where the i-th probe colour was detected in the image.
%   `probe_regions_bw` would normally be produced by
%   'detectWithinBoundingBox()', as its `probe_regions_bw_filtered` output
%   argument.
%
% model_pca_space -- Bilateral model of detected probe edge endpoints (PCA space)
%   The `model` output argument of 'bilateralModel()', when
%   'bilateralModel()' is called on the probe band edge endpoints detected
%   from the regions in `probe_regions_bw`.
%
% model -- Bilateral model of detected probe edge endpoints (image space)
%   The `model_px` output argument of 'bilateralModel()', when
%   'bilateralModel()' is called on the probe band edge endpoints detected
%   from the regions in `probe_regions_bw`.
%
% model_transform -- Mapping from PCA space to image space
%   The `transform` output argument of 'bilateralModel()', when
%   'bilateralModel()' is called on the probe band edge endpoints detected
%   from the regions in `probe_regions_bw`.
%
% params -- Fixed parameters
%   Parameters that should be stable across a variety of input images and
%   probe models. `params` is a structure containing the following fields:
%   - color_sum_outlier_threshold: Number of standard deviations from
%       the average score for a colour label below which a colour label is
%       deemed to have insufficient support. The region is assigned a label
%       of 'no colour' (0).
%   - color_dominance_threshold: If the colour label with the highest score
%       for a given region has a score that is less than
%       `color_dominance_threshold` times the score of the label with the
%       second-highest score, then the region is assigned a label of 'any
%       colour' (-1).
%
% verbose -- Debugging flags
%   If recognized fields of `verbose` are true, corresponding graphical
%   output will be generated for debugging purposes.
%
%   All debugging flags default to false if `verbose` is not passed.
%
% ## Output Arguments
%
% colors_left -- Labels for colours to the left of detected band edges
%   A column vector containing the labels assigned to the colours to the
%   left of the detected edges between the coloured bands of the probe.
%   "left" means smaller coordinates along the first PCA axis defined by
%   `model_pca_space` than the coordinates of the detected edges.
%
% colors_right -- Labels for colours to the right of detected band edges
%   Similar to `colors_left`, but for colours to the right of the detected
%   edges between the coloured bands of the probe.
%
% ## Algorithm
%
% The probe has a series of coloured bands, belonging to two or more colour
% classes, as labelled by the user. Edges between bands therefore lie
% between regions with different colour labels. It may seem trivial to
% determine the labels on either side of a detected edge, given that edges
% are detected by finding narrow spaces between detected regions with
% different colour labels. However, some of the detected regions are false
% positives. Other detected regions may poorly overlap with their ground
% truth shapes in the image.
%
% This function attempts to robustly determine the appropriate colour label
% between detected edges, by considering the total area between the edges
% occupied by each colour. The score for a given colour is a weighted count
% of the pixels assigned to that colour between two detected edges. Weights
% are presently distances (measured parallel to the estimated axis of the
% probe) from the edge for which the colour label is being computed. Only
% pixels lying within the bounding box defined by the edge endpoints are
% counted.
%
% Following the weighted counting operation, colour label scores are
% processed according to the fields of `params`, as described above.
%
% ## Notes
% - This function computes binary images for the probe colours that are
%   clipped to the bounding quadrilaterals defined by the image edges and
%   by the detected endpoints of the edges between probe colour bands.
%   These binary images might be useful for other tasks, and could be
%   assigned to an output argument in the future.
%
% See also detectWithinBoundingBox, bilateralModel

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created March 28, 2017

%% Parse input arguments

nargoutchk(2,2);
narginchk(8,9);

if ~isempty(varargin)
    verbose = varargin{1};
    display_final_clipped_regions_colored = verbose.display_final_clipped_regions_colored;
else
    display_final_clipped_regions_colored = false;
end

image_width = size(I, 2);
image_height = size(I, 1);

n_detected_band_edges = size(model_pca_space.above, 1);

%% Label the colours between detected probe band edges

% Clip the detected colour regions to the area between the probe edge
% endpoints. Extrapolate the clipping area to the tips of the probe by extending the
% line segments formed by the first and last pairs of band edges.

% Interior region
polygons_points = zeros((n_detected_band_edges + 2) * 2, 2);
polygons_points(2:(n_detected_band_edges + 1), :) = model.above;
polygons_points((n_detected_band_edges + 4):(end-1), :) = flipud(model.below);

% Define line segments used to extrapolate to the probe tips
line_endpoints = [
    model.above(1,:), model.above(2,:);
    model.above(end,:), model.above(end-1,:);
    model.below(end,:), model.below(end-1,:);
    model.below(1,:), model.below(2,:)
    ];

% Find intersections with the four image borders
% Top, right, bottom, left
border_coordinates = [1; image_width; image_height; 1];
border_coordinate_indices = [2; 1; 2; 1];

line_endpoints_rep = repmat(line_endpoints, 4, 1);
border_coordinates_rep = repelem(border_coordinates, 4);
border_coordinate_indices = repelem(border_coordinate_indices, 4);
line_indices = repmat((1:size(line_endpoints, 1))', 4, 1);

border_points = pointsOnLinesByCoordinates(...
    line_endpoints_rep(:, 1:2), line_endpoints_rep(:, 3:4),...
    border_coordinates_rep, border_coordinate_indices);

% Determine which intersection points to use
border_points_filter = all(isfinite(border_points), 2);
% Select intersection points in the appropriate directions
border_points_filter = border_points_filter &...
    (dot(border_points - line_endpoints_rep(:, 1:2),...
     line_endpoints_rep(:, 1:2) - line_endpoints_rep(:, 3:4), 2) > 0);
% Select intersection points within the image
border_points_filter = border_points_filter &...
    border_points(:, 1) > 0 & border_points(:, 1) <= image_width &...
    border_points(:, 2) > 0 & border_points(:, 2) <= image_height;
line_indices = line_indices(border_points_filter);
border_points = border_points(border_points_filter, :);
border_points(line_indices, :) = border_points;

% If lines crossed, replace border points by line intersections
if dot(border_points(1, :) - border_points(4, :), line_endpoints(1, 1:2) - line_endpoints(4, 1:2), 2) < 0
    intersection_point = cross(...
            cross([border_points(1, :), 1], [line_endpoints(1, 1:2), 1]),...
            cross([border_points(4, :), 1], [line_endpoints(4, 1:2), 1])...
        );
    border_points([1 4], :) = repmat(intersection_point(1:2) ./ intersection_point(3), 2, 1);
end
if dot(border_points(2, :) - border_points(3, :), line_endpoints(2, 1:2) - line_endpoints(3, 1:2), 2) < 0
    intersection_point = cross(...
            cross([border_points(2, :), 1], [line_endpoints(2, 1:2), 1]),...
            cross([border_points(3, :), 1], [line_endpoints(3, 1:2), 1])...
        );
    border_points([2 3], :) = repmat(intersection_point(1:2) ./ intersection_point(3), 2, 1);
end

polygons_points([1, n_detected_band_edges + 2, n_detected_band_edges + 3, end], :) = border_points;

% Perform the clipping
clipping_mask = roipoly(I, polygons_points(:, 1), polygons_points(:, 2));
probe_regions_bw_final_clipped = probe_regions_bw;

n_colors = size(probe_regions_bw, 3);
for i = 1:n_colors
    probe_regions_bw_final_clipped(:, :, i) = probe_regions_bw_final_clipped(:, :, i) & clipping_mask;
end

if display_final_clipped_regions_colored
    probe_regions_bw_final_clipped_display = zeros(image_height, image_width, 3);
    for i = 1:n_colors
        [ ~, peak_hue_index ] = max(probe_color_distributions(:, i));
        peak_hue = (peak_hue_index - 1) * probe_color_distribution_increment;
        peak_rgb = hsv2rgb([peak_hue, 1, 1]);
        probe_regions_bw_final_clipped_filtered_i = probe_regions_bw_final_clipped(:, :, i);
        probe_regions_bw_final_clipped_display = probe_regions_bw_final_clipped_display +...
            cat(3,...
                    peak_rgb(1) * probe_regions_bw_final_clipped_filtered_i,...
                    peak_rgb(2) * probe_regions_bw_final_clipped_filtered_i,...
                    peak_rgb(3) * probe_regions_bw_final_clipped_filtered_i...
                );
    end
    figure
    imshow(probe_regions_bw_final_clipped_display);
    title('Final detected probe regions, clipped to detected probe edges')
end

% Find the first components of the PCA-space coordinates of the probe colour regions
probe_regions_pca_space_x = cell(n_colors, 1);
for i = 1:n_colors
    probe_regions_bw_final_clipped_i = probe_regions_bw_final_clipped(:, :, i);
    pixel_indices_i = find(probe_regions_bw_final_clipped_i);
    [y,x] = ind2sub([image_height, image_width],pixel_indices_i);
    probe_regions_pca_space_xi = [x, y, ones(length(x),1)];
    probe_regions_pca_space_xi = probe_regions_pca_space_xi / model_transform;
    % Normalization (dividing by the homogenous coordinate) is not
    % necessary, as `model_transform` is an affine transformation.
    probe_regions_pca_space_x{i} = probe_regions_pca_space_xi(:, 1);
end

% Iterate over each probe band to find the colours dominating at either end
probe_color_scores_left = zeros(n_detected_band_edges, n_colors);
probe_color_scores_right = zeros(n_detected_band_edges, n_colors);
for i = 1:n_detected_band_edges
    if i > 1
        min_coord = mean(...
            [model_pca_space.above(i-1,1),...
            model_pca_space.below(i-1,1)]...
            );
    else
        min_coord = -Inf;
    end
    if i < n_detected_band_edges
        max_coord = mean(...
            [model_pca_space.above(i+1,1),...
            model_pca_space.below(i+1,1)]...
            );
    else
        max_coord = Inf;
    end
    center_coord = mean(...
            [model_pca_space.above(i,1),...
            model_pca_space.below(i,1)]...
            );
    for j = 1:n_colors
        probe_regions_pca_space_xi = probe_regions_pca_space_x{j};
        left_x = probe_regions_pca_space_xi > min_coord &...
            probe_regions_pca_space_xi <= center_coord;
        left_x = probe_regions_pca_space_xi(left_x);
        left_x = repmat(center_coord, length(left_x), 1) - left_x;
        probe_color_scores_left(i,j) = sum(left_x);
        right_x = probe_regions_pca_space_xi >= center_coord &...
            probe_regions_pca_space_xi < max_coord;
        right_x = probe_regions_pca_space_xi(right_x);
        right_x = right_x - repmat(center_coord, length(right_x), 1);
        probe_color_scores_right(i,j) = sum(right_x);
    end
end

% Pick the dominant colour, provided it dominates by a sufficient factor
for i = 1:2
    if i == 1
        scores = probe_color_scores_left;
    else
        scores = probe_color_scores_right;
    end
    [scores_max, scores_maxind] = max(scores, [], 2);
    
    % Filter out outlier small scores
    mu_scores_max = mean(scores_max);
    sigma_scores_max = std(scores_max);
    mu_scores_max_rep = repmat(mu_scores_max,n_detected_band_edges,1);
    sigma_scores_max_rep = repmat(sigma_scores_max,n_detected_band_edges,1);
    scores_max_filter = (mu_scores_max_rep - scores_max) > params.color_sum_outlier_threshold * sigma_scores_max_rep;
    scores_max(scores_max_filter) = 0;

    scores_nomax = scores;
    scores_nomax(...
        sub2ind(size(scores), (1:n_detected_band_edges)', scores_maxind)...
        ) = 0;
    scores_max2 = max(scores_nomax, [], 2);
    colors = scores_max ./ scores_max2;
    colors_filter = colors >= params.color_dominance_threshold;
    colors(~colors_filter) = -1; % Any colour
    colors(scores_max == 0) = 0; % No colour
    colors(colors_filter) = scores_maxind(colors_filter);
    if i == 1
        colors_left = colors;
    else
        colors_right = colors;
    end
end

end