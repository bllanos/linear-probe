function [ colors_left, colors_right ] = labelDetectedColors(...
    I, probe_color_distributions,...
    probe_regions_bw,...
    model_pca_space, model, model_transform, axis, varargin...
)
%LABELDETECTEDCOLORS Robustly associate detected probe regions with colour classes
%
% ## Syntax
% [ colors_left, colors_right ] = labelDetectedColors(...
%     I, probe_color_distributions,...
%     probe_regions_bw,...
%     model_pca_space, model, model_transform, axis [, verbose]...
% )
%
% ## Description
% [ colors_left, colors_right ] = labelDetectedColors(...
%     I, probe_color_distributions,...
%     probe_regions_bw,...
%     model_pca_space, model, model_transform, axis [, verbose]...
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
%   Discretized density estimators of image hue values corresponding to the
%   different coloured bands on the probe, in the same order (starting from
%   the active tip of the probe). The i-th column of this 2D array stores
%   the estimator for the i-th colour class of probe segments. (There are
%   `n_colors` columns in total.)
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
% axis -- Estimated probe axis in image
%   A 3-vector containing the homogenous coordinates of a line in the image
%   approximating the probe axis.
%
% verbose -- Debugging flag
%   If true, graphical output will be generated for debugging purposes.
%
%   Defaults to false if not passed.
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
% different colour labels. However, detected regions may poorly overlap
% with their ground truth shapes in the image, and there may be spurious,
% small detected regions.
%
% Instead of assigning labels based on the closest detected colour regions
% to each detected edge, this function assigns labels based on the closest
% centroids of detected colour regions to each detected edge. Centroids are
% computed after clipping the detected coloured regions to the
% quadrilaterals defined by the endpoints of adjacent detected edges.
% Assigning colour labels based on distances to centroids provides
% robustness to slight biases of detected edge positions towards the
% interiors of detected colour regions.
%
% Centroids are first filtered to those whose regions (as they appear after
% clipping) cross the estimated axis of the probe. This step eliminates
% many small, spurious detected regions.
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

% Parse input arguments

nargoutchk(2,2);
narginchk(7,8);

if ~isempty(varargin)
    verbose = varargin{1};
else
    verbose = false;
end

image_width = size(I, 2);
image_height = size(I, 1);

n_detected_band_edges = size(model_pca_space.above, 1);

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
     line_endpoints_rep(:, 1:2) - line_endpoints_rep(:, 3:4), 2) >= 0);
% Select intersection points within the image
border_points_filter = border_points_filter &...
    border_points(:, 1) >= 1 & border_points(:, 1) <= image_width &...
    border_points(:, 2) >= 1 & border_points(:, 2) <= image_height;
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
probe_regions_bw_clipped = probe_regions_bw;

n_colors = size(probe_regions_bw, 3);
for i = 1:n_colors
    probe_regions_bw_clipped(:, :, i) = probe_regions_bw_clipped(:, :, i) & clipping_mask;
end

if verbose
    probe_regions_bw_clipped_display = zeros(image_height, image_width, 3);
    probe_color_distribution_increment = hueSamplingParams( probe_color_distributions(:, 1) );
    for i = 1:n_colors
        [ ~, peak_hue_index ] = max(probe_color_distributions(:, i));
        peak_hue = (peak_hue_index - 1) * probe_color_distribution_increment;
        peak_rgb = hsv2rgb([peak_hue, 1, 1]);
        probe_regions_bw_clipped_filtered_i = probe_regions_bw_clipped(:, :, i);
        probe_regions_bw_clipped_display = probe_regions_bw_clipped_display +...
            cat(3,...
                    peak_rgb(1) * probe_regions_bw_clipped_filtered_i,...
                    peak_rgb(2) * probe_regions_bw_clipped_filtered_i,...
                    peak_rgb(3) * probe_regions_bw_clipped_filtered_i...
                );
    end
    figure
    imshow(probe_regions_bw_clipped_display);
    title('Detected probe regions, clipped to detected probe edges')
end

% Extract region centroids

regions = struct('Connectivity', cell(n_colors, 1), 'ImageSize', cell(n_colors, 1),...
    'NumObjects', cell(n_colors, 1), 'PixelIdxList', cell(n_colors, 1));

for i = 1:n_colors
    regions(i) = bwconncomp(probe_regions_bw_clipped(:, :, i));
end

centroids = cell(n_colors, 1);
for i = 1:n_colors
    s = regionprops(regions(i), 'centroid');
    centroids{i} = cat(1, s.Centroid);
end
centroids_all = cell2mat(centroids);
n_regions_all = size(centroids_all, 1);

centroids_color_index = zeros(n_regions_all, 1);
offset = 0;
for i = 1:n_colors
    centroids_color_index((1+offset):(offset + size(centroids{i}, 1))) = i;
    offset = offset + size(centroids{i}, 1);
end

% Filter out regions which do not cross the estimated probe axis

% Find region pixel coordinates
image_points = cell(n_regions_all, 1);
offset = 1;
for i = 1:n_colors
    px_indices = regions(i).PixelIdxList;
    for j = 1:regions(i).NumObjects
        [row, col] = ind2sub([image_height, image_width], px_indices{j});
        image_points{offset} = [col, row];
        offset = offset + 1;
    end
end

% Filter to regions which cross the estimated probe axis
regions_filter = isOnBothSidesOfLine( image_points, axis );
centroids_filtered = centroids_all(regions_filter, :);
centroids_color_index = centroids_color_index(regions_filter);
n_regions_filtered = length(centroids_color_index);

% Label the colours between detected probe band edges

% Find the first components of the PCA-space coordinates of the centriods
centroids_pca_space_x = [
        centroids_filtered(:, 1),...
        centroids_filtered(:, 2),...
        ones(n_regions_filtered,1)
        ] / model_transform;
% Normalization (dividing by the homogenous coordinate) is not
% necessary, as `model_transform` is an affine transformation.
centroids_pca_space_x = centroids_pca_space_x(:, 1);

% Iterate over each probe band to find the colour of the closest centroid
colors_left = zeros(n_detected_band_edges, 1);
colors_right = zeros(n_detected_band_edges, 1);
for i = 1:n_detected_band_edges
    center_coord = mean(...
            [model_pca_space.above(i,1),...
            model_pca_space.below(i,1)]...
            );
    
    distances = centroids_pca_space_x - center_coord;
    distances_left = distances;
    distances_left(distances > 0) = -Inf;
    [val, index] = max(distances_left);
    if isfinite(val)
        colors_left(i) = centroids_color_index(index);
    end
    
    distances_right = distances;
    distances_right(distances < 0) = Inf;
    [val, index] = min(distances_right);
    if isfinite(val)
        colors_right(i) = centroids_color_index(index);
    end
end

if verbose    
    figure
    imshow(probe_regions_bw_clipped_display);
    hold on
    scatter(centroids_all(~regions_filter, 1), centroids_all(~regions_filter, 2), 'ro');
    scatter(centroids_filtered(:, 1), centroids_filtered(:, 2), 'go');
    hold off
    title('Centroids used for colour labelling (green) and rejected centriods (red)')
end

end