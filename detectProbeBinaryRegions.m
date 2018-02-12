function [ regions_filtered, bw_filtered ] = detectProbeBinaryRegions( regions, bw, color_pairs, radius_adj, axis_distance_outlier_threshold, varargin )
% DETECTPROBEBINARYREGIONS  Filter binary regions to those which could belong to the probe
%
% ## Syntax
% regions_filtered = detectProbeBinaryRegions(...
%   regions, bw, color_pairs, radius_adj,...
%   axis_distance_outlier_threshold [, verbose]...
% )
% [ regions_filtered, bw_filtered ] = detectProbeBinaryRegions(...
%   regions, bw, color_pairs, radius_adj,...
%   axis_distance_outlier_threshold [, verbose]...
% )
%
% ## Description
% regions_filtered = detectProbeBinaryRegions(...
%   regions, bw, color_pairs, radius_adj,...
%   axis_distance_outlier_threshold [, verbose]...
% )
%   Returns binary regions that satisfy basic probe geometry constraints.
%
% [ regions_filtered, bw_filtered ] = detectProbeBinaryRegions(...
%   regions, bw, color_pairs, radius_adj,...
%   axis_distance_outlier_threshold [, verbose]...
% )
%   Additionally returns binary images corresponding to the binary regions.
%
% ## Input Arguments
%
% regions -- Initial binary regions
%   A structure vector of length n, where the i-th element contains the
%   candidate regions corresponding to the i-th colour. Each element is of
%   the form of the output argument of `bwconncomp`.
%
% bw -- Initial binary images
%   An image_height x image_width x n logical array, where `bw(:,:,i)` is a
%   binary image containing the regions tentatively belonging to the i-th
%   colour class of the probe, contained in `regions(i)`.
%
% color_pairs -- Colour adjacency relationships
%   A k x 2 array, where each row contains the indices of two colour
%   classes of the probe. A row appears in `color_pairs` if there are two
%   adjacent bands on the probe with the corresponding colours to the row.
%   (Bands can be in any order.)
%
%   `color_pairs` should not contain any duplicate pairs.
%
% radius_adj -- Radius for region proximity test
%   Threshold above which regions are rejected based on their smallest
%   distance to any region belonging to a different colour class.
%   Specifically, the adjacent regions have the colour classes with which
%   the colour class of the given region is paired in `colour_pairs`.
%
%   This test enforces the requirement that the coloured bands of the probe
%   are adjacent. Note that bands next to occluding objects may be rejected
%   by this test, but this is acceptable, because no usable endpoints can
%   be extracted from such bands for the purposes of probe pose estimation.
%
% axis_distance_outlier_threshold - Outlier threshold in standard deviations
%   The number of standard deviations from the mean perpendicular distance
%   to the estimated axis of the probe beyond which a region is rejected
%   based on the position of its centroid.
%
%   This test enforces the requirement that the coloured bands of the probe
%   form a line.
%
% verbose -- Debugging flag
%   If true, graphical output will be generated for debugging purposes.
%
%   Defaults to false if not passed.
%
% ## Output Arguments
%
% regions_filtered -- Filtered binary regions
%   A structure vector of length n, where the i-th element contains the
%   regions corresponding to the i-th colour that pass all tests. Each
%   element is of the form of the output argument of `bwconncomp`.
%
% bw_filtered - Filtered binary images
%   An image_height x image_width x n logical array, where
%   `bw_filtered(:,:,i)` is the binary image consisting of the connected
%   components in `regions_filtered(i)`.
%
% ## References
% Section 4.7 (on RANSAC) of Hartley, Richard, and Andrew Zisserman. Multiple
%   View Geometry In Computer Vision. Cambridge, UK: Cambridge University
%   Press, 2003. eBook Academic Collection (EBSCOhost). Web. 15 August
%   2016.
%
% See also bwconncomp, extractBinaryRegions

% Bernard Llanos
% Spring 2016 research assistantship supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created August 15, 2016

nargoutchk(1, 2);
narginchk(5, 6);

if ~isempty(varargin)
    verbose = varargin{1};
else
    verbose = false;
end

image_height = size(bw, 1);
image_width = size(bw, 2);
n = size(bw, 3);

regions_filtered = struct('Connectivity', cell(n, 1), 'ImageSize', cell(n, 1),...
    'NumObjects', cell(n, 1), 'PixelIdxList', cell(n, 1));

bw_filtered_perim = false(image_height, image_width, n);

% Filter regions to those within `radius_adj` of regions for a possible
% adjacent colour
for i = 1:n
    adjacent_colors = [
            color_pairs(color_pairs(:, 1) == i, 2);
            color_pairs(color_pairs(:, 2) == i, 1)
        ];
    % Binary image containing all regions for all possible adjacent colours
    bw_adj_i = any(bw(:, :, adjacent_colors), 3);
    bw_distance = bwdist(bw_adj_i);
    bw_distance = (bw_distance <= radius_adj);
    if verbose
        figure
        imshow(bw(:, :, i) & bw_distance);
        title(sprintf('Binary image %d AND-ed with thresholded distances to regions for adjacent colours', i))
    end
    bw_filtered_perim_i = false(image_height, image_width);
    regions_i = regions(i);
    region_pixel_lists_i = regions_i.PixelIdxList;
    distance_filter = false(regions_i.NumObjects, 1);
    for j = 1:regions_i.NumObjects
        if any(bw_distance(region_pixel_lists_i{j}))
            distance_filter(j) = true;
            bw_filtered_perim_i(region_pixel_lists_i{j}) = true;
        end
    end
    regions_i.PixelIdxList = region_pixel_lists_i(distance_filter);
    regions_i.NumObjects = length(regions_i.PixelIdxList);
    bw_filtered_perim(:, :, i) = bwperim(bw_filtered_perim_i);
    regions_filtered(i) = regions_i;
end

% Filter regions to those close to the estimated axis of the probe

% Find region centroids
centroids = cell(n, 1);
for i = 1:n
    s = regionprops(regions_filtered(i), 'centroid');
    centroids{i} = cat(1, s.Centroid);
end
centroids = cell2mat(centroids);
n_regions_all = size(centroids, 1);

if n_regions_all < 2
    error('PROBE:InsufficientRegionsDetected', 'Less than two regions within %g pixels of an appropriate adjacent colour were detected.', radius_adj)
end

% Find region pixel coordinates
image_points = cell(n_regions_all, 1);
offset = 1;
for i = 1:n
    px_indices = regions_filtered(i).PixelIdxList;
    bw_filtered_perim_i = bw_filtered_perim(:, :, i);
    for j = 1:regions_filtered(i).NumObjects
        % Only use points on the perimeter
        perim_filter = bw_filtered_perim_i(px_indices{j});
        [row, col] = ind2sub([image_height, image_width], px_indices{j}(perim_filter));
        image_points{offset} = [col, row];
        offset = offset + 1;
    end
end

% RANSAC-based estimation of probe axis
trial = 1;
n_trials = nchoosek(n_regions_all, 2); % Number of outliers is unknown - This is the worst case
sample_indices = nchoosek(1:n_regions_all, 2);
sample_indices = sample_indices(randperm(n_trials), :);
n_inliers_max = 0;
inliers_count_threshold = n_regions_all; % Number of inliers is unknown - This is the worst case
p = 0.99; % Require this probability that at least one of `n_trials` samples contains only inliers
while trial <= n_trials
    sample = centroids(sample_indices(trial, :), :);
    
    % Rather than define outliers based on their distance to the estimated
    % axis (as in past versions of this function), define outliers as
    % regions which contain pixels on only one side of the estimated axis.
    inliers_filter_candidate = isOnBothSidesOfLine( image_points, sample );
    
    n_inliers = sum(inliers_filter_candidate);
    if n_inliers > n_inliers_max
        successful_sample = sample;
        n_inliers_max = n_inliers;
        inliers_filter = inliers_filter_candidate;
        n_trials = min(n_trials, log(1 - p) / log(1 - (n_inliers_max / n_regions_all) ^ 2));
        if n_inliers >= inliers_count_threshold
            break;
        end
    end
    trial = trial + 1;
end
n_inliers = n_inliers_max;

% Finalize set of inliers
% This time, use distances from region centroids to the estimated probe
% axis, which is a less restrictive outlier detection method than the one
% used during RANSAC iteration. The previous outlier detection method is
% too harsh, given the uncertainty in the estimated probe axis.
diff_n_inliers = 1;
while n_inliers > 1 && diff_n_inliers > 0
    [coeff_pca, ~, ~, ~, ~, mu_pca] = pca(centroids(inliers_filter, :));
    centroids_pca = (centroids - repmat(mu_pca, n_regions_all, 1)) / (coeff_pca.');
    if size(centroids_pca, 2) == 1 % Centroids span a 1D space
        distances = zeros(length(centroids_pca), 1);
    else
        distances = abs(centroids_pca(:, 2));
    end
    mu = mean(distances(inliers_filter));
    sigma = std(distances(inliers_filter));
    mu_rep = repmat(mu,n_regions_all,1);
    sigma_rep = repmat(sigma,n_regions_all,1);
    inliers_filter = abs(distances - mu_rep) <= axis_distance_outlier_threshold * sigma_rep;
    n_inliers_new = sum(inliers_filter);
    diff_n_inliers = abs(n_inliers_new - n_inliers);
    n_inliers = n_inliers_new;
end

% Filter regions to inlier regions
offset = 1;
for i = 1:n
    regions_filtered(i).PixelIdxList = regions_filtered(i).PixelIdxList(...
        inliers_filter(offset:(offset + regions_filtered(i).NumObjects - 1))...
        );
    offset = offset + regions_filtered(i).NumObjects;
    regions_filtered(i).NumObjects = length(regions_filtered(i).PixelIdxList);
end

if nargout > 1
    bw_filtered = false(image_height, image_width, n);
    for i = 1:n
        regions_filtered_i = regions_filtered(i);
        region_filtered_pixel_lists_i = regions_filtered_i.PixelIdxList;
        bw_filtered_i = false(image_height, image_width);
        for j = 1:regions_filtered_i.NumObjects
            bw_filtered_i(region_filtered_pixel_lists_i{j}) = true;
        end
        bw_filtered(:, :, i) = bw_filtered_i;
        if verbose
            figure
            imshow(bw_filtered_i);
            title(sprintf('Final binary image %d', i))
            hold on
            line(successful_sample(:, 1), successful_sample(:, 2))
            hold off
        end
    end
end

end

