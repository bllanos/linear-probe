function [ endpoints ] = detectProbeEdgeEndpoints( bw, edge_width, sigma_theta, threshold, varargin )
% DETECTPROBEEDGEENDPOINTS  Estimate where edges between probe bands intersect the outer edges of the probe
%
% ## Syntax
% endpoints = detectProbeEdgeEndpoints( bw, edge_width, sigma_theta, threshold [, verbose] )
%
% ## Description
% endpoints = detectProbeEdgeEndpoints( bw, edge_width, sigma_theta, threshold [, verbose] )
%   Returns the pixel coordinates of the estimated locations of the visible
%   endpoints of the input edges between probe colour regions.
%
% ## Input Arguments
%
% bw -- Probe colour region junctions
%   An image_height x image_width binary image containing approximate edges
%   between the coloured bands of the probe.
%
% edge_width -- Band edge width
%   The pixel width of the edges in `bw`, perhaps better defined as the
%   search distance used to obtain the edges. This parameter is used to set
%   the falloff of the radial component of the filter used to refine the
%   edges in `bw` to join broken edges and eliminate portions of edges that
%   protrude along the outside the probe. It is also used to eliminate
%   short edges.
%
% sigma_theta -- Angular filter width
%   The standard deviation of the angular component of the filter used to
%   refine the edges in `bw`. A smaller value results in greater
%   suppression of edges that deviate from the orientation of the minor
%   axis of the probe.
%
% threshold -- Edge refinement threshold
%   The threshold used to binarize the result of filtering the input edge
%   image, `bw`, with a filter that emphasizes portions of edges aligned
%   with the minor axis of the probe.
%
% verbose -- Debugging flag
%   If true, graphical output will be generated for debugging purposes.
%
%   Defaults to false if not passed.
%
% ## Output Arguments
%
% endpoints -- Estimated edge endpoints
%   A k x 2 array, where the columns store the x and y pixel coordinates of
%   the final set of k edge endpoints. This is the set of points estimated
%   to be located at the endpoints of the edges between the coloured bands
%   of the probe - The points where the edges vanish at the outer borders
%   of the probe in the image.
%
% ## Algorithm
% 
% The input image, `bw`, is assumed to contain binary regions which are
% fairly noisy approximations to the edges between the bands of the probe.
% Specifically, regions may spill outwards along the exterior of the probe
% and may be broken due to strong highlights and reflections, but are
% assumed to at least contain the true points where the edges meet the
% outer borders of the probe. (This function cannot compensate for edges
% that have been detected shorter than they actually are. In general, there
% is no way to correct for this without using information about the
% physical widths of the junctions between bands, which can only be used
% after the detected edges have been matched to probe band junctions.)
%
% The process used to estimate where band edges intersect the probe's
% borders is as follows:
% - The probe's orientation is determined by running Principal Components
%   Analysis on the regions in the input edge image.
%
% - The input edge image is filtered with a filter kernel that suppresses
%   portions of edges misaligned with the probe's minor axis. The
%   `edge_width` parameter is the standard deviation of the radial Gaussian
%   component of the filter, whereas the `sigma_theta` parameter is the
%   standard deviation of the angular Gaussian component of the filter.
%
% - The filtered image is binarized with `threshold` and intersected with
%   the original edge image. The resulting binary image should contain
%   fewer portions of edges that protrude along the outer borders of the
%   probe.
%
% - Candidate edge endpoints are extracted by finding the extreme points in
%   each binary region along the PCA minor axis of the probe.
%
% - As there may still be spurious edge fragments lying along the outer
%   borders of the probe, edge endpoints are rejected if their separation
%   is less than `edge_width`. Note that this may incorrectly reject small
%   portions of edges broken by highlights or reflections on the probe.
%   This is not critical as long as these edge portions are not located
%   where the true edges intersect the borders of the probe.
%
% - Candidate edge endpoints are removed if they are not on the expected
%   side of the PCA major axis. Specifically, the endpoint with the largest
%   second coordinate in the space of the PCA axes should be above the
%   major axis, whereas the other endpoint should be below the major axis.
%   This step will incorrectly reject endpoints if the PCA major axis is
%   located outside the true borders of the probe.
%
% - At this stage, binary regions which correctly span the width of the
%   probe each provide two final endpoints. However, binary regions forming
%   broken edges still have spurious endpoints located inside the probe's
%   borders (although fewer than before the preceding steps). It may be
%   possible to remove these spurious endpoints by clustering endpoints
%   into groups approximately aligned with the probe's minor axis, but this
%   will incorrectly group endpoints from the two sides of very thin bands.
%   Instead, the function attempts to repair broken edges by operating on
%   the binary image obtained by filtering above (before it was intersected
%   with the original edge image). For each set of candidate edge points
%   located in the same binary region, only those with the most extreme
%   second coordinates in the PCA axes space are retained.
%
% See also pca, bilateralModel

% Bernard Llanos
% Spring 2016 research assistantship supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created August 22, 2016

    function [ H ] = directionalKernel(s_distance, s_theta, mu_theta)
        % Returns a filtering kernel which responds to pixels along lines
        % of high intensity with orientations close to `mu_theta`. Filter
        % weights are obtained from a Gaussian function of distance from
        % the center of the kernel and a Gaussian function of the angle of
        % the neighbouring pixel relative to `mu_theta`.
        %
        % Angles are expressed in radians, in the range -pi/2 to pi/2 from
        % the positive x-axis.
        
        % Define the size of the filter so that the weight at its outer
        % edge is at most `p`.
        p = 0.025 / normpdf(0, 0, s_theta);
        radial_limit = norminv([p p], 0, s_distance);
        radial_limit = ceil(abs(radial_limit(1)));
        H = zeros(radial_limit, radial_limit);
        center = floor((size(H) + 1) / 2);
        for i_inner = 1:radial_limit
            for j_inner = 1:radial_limit
                dx = j_inner - center(2);
                dy = i_inner - center(1);
                r = sqrt(dx ^ 2 + dy ^ 2);
                theta = atan(dy / dx);
                theta_deviation = abs(theta - mu_theta);
                theta_deviation = min(theta_deviation, pi - theta_deviation);
                H(i_inner, j_inner) = normpdf(r, 0, s_distance) * normpdf(theta_deviation, 0, s_theta);
            end
        end
        H(center(1), center(2)) = normpdf(0, 0, s_distance) * normpdf(0, 0, s_theta);
        % Normalize weights
        H = H ./ sum(sum(H));
    end

nargoutchk(1, 1);
narginchk(4, 5);

if ~isempty(varargin)
    verbose = varargin{1};
else
    verbose = false;
end

image_height = size(bw, 1);
image_width = size(bw, 2);
image_size = [image_height, image_width];

% Construct a filtering kernel to suppress lines parallel to the probe's
% major axis
[edge_points_y, edge_points_x] = find(bw);
coeff_pca = pca([edge_points_x, edge_points_y]);
pca_angle = atan(coeff_pca(2, 2) / coeff_pca(1, 2));
kernel = directionalKernel(...
    edge_width, sigma_theta, pca_angle);

if verbose
    figure;
    surf(kernel);
    title('Filter kernel for suppressing lines parallel to the probe edges')
    xlabel('X')
    ylabel('Y')
    zlabel('Weight')
end

% Filter the image and threshold it to obtain refined edges
bw_double = im2double(bw);
bw_filtered = imfilter(bw_double, kernel);
if verbose
    figure;
    bw_comparison = cat(3, bw_double, bw_filtered, zeros(image_height, image_width));
    imshow(bw_comparison);
    title('Filtered edge image (green) and original (red)')
end
bw_filtered = imbinarize(bw_filtered, threshold);
bw_new = bw_filtered & bw;

% Find edges and their endpoints along the second PCA axis
[edge_points_y, edge_points_x] = find(bw_new);
[~,score] = pca([edge_points_x, edge_points_y]);
edge_points_linear = sub2ind(image_size, edge_points_y, edge_points_x);
score_image = zeros(image_height, image_width, 2);
score_image(edge_points_linear) = score(:, 1);
score_image(edge_points_linear + image_height * image_width) = score(:, 2);
edge_points_by_edge = bwconncomp(bw_new);
edge_endpoints_pca = zeros(edge_points_by_edge.NumObjects * 2, 2);
edge_endpoints = zeros(edge_points_by_edge.NumObjects * 2, 2);
k = 1;
for i = 1:edge_points_by_edge.NumObjects
    edge_points_i = edge_points_by_edge.PixelIdxList{i};
    edge_points_pca_2 = score_image(edge_points_i + image_height * image_width);
    [~, min_ind] = min(edge_points_pca_2);
    [~, max_ind] = max(edge_points_pca_2);
    edge_endpoints_pca(k, :) = [
            score_image(edge_points_i(min_ind)),...
            edge_points_pca_2(min_ind)
        ];
    edge_endpoints_pca(k + 1, :) = [
            score_image(edge_points_i(max_ind)),...
            edge_points_pca_2(max_ind)
        ];
    [ edge_endpoints(k, 2), edge_endpoints(k, 1) ] =...
        ind2sub(image_size, edge_points_i(min_ind));
    [ edge_endpoints(k + 1, 2), edge_endpoints(k + 1, 1) ] =...
        ind2sub(image_size, edge_points_i(max_ind));
    k = k + 2;
end

% Remove endpoints on the same edges having a separation on the order of
% the edge width
edge_endpoints_separation = diff(edge_endpoints, 1);
edge_endpoints_separation = edge_endpoints_separation(1:2:end, :);
edge_endpoints_separation = dot(edge_endpoints_separation, edge_endpoints_separation, 2);
edge_endpoints_filter = edge_endpoints_separation > (edge_width  ^ 2);
edge_endpoints_filter = repelem(edge_endpoints_filter, 2);
edge_endpoints_filtered = edge_endpoints(edge_endpoints_filter, :);
edge_endpoints_pca_filtered = edge_endpoints_pca(edge_endpoints_filter, :);

% Remove endpoints on the same edges that are not separated by the PCA
% major axis of the probe
edge_endpoints_filter = false(size(edge_endpoints_filtered, 1), 1);
edge_endpoints_filter(1:2:end) = edge_endpoints_pca_filtered(1:2:end, 2) < 0;
edge_endpoints_filter(2:2:end) = edge_endpoints_pca_filtered(2:2:end, 2) > 0;
edge_endpoints_filtered = edge_endpoints_filtered(edge_endpoints_filter, :);
edge_endpoints_pca_filtered = edge_endpoints_pca_filtered(edge_endpoints_filter, :);

% Attempt to join broken edges
%
% Find edge points in the same binary regions in the filtered image
[ bw_filtered_labelled, n_labels ] = bwlabel(bw_filtered);
edge_endpoints_labels = bw_filtered_labelled(...
        sub2ind(image_size, edge_endpoints_filtered(:, 2), edge_endpoints_filtered(:, 1))...
    );
% Select only the minimum and maximum points along the PCA minor axis for
% each binary region
edge_endpoints_filter = false(size(edge_endpoints_filtered, 1), 1);
for i = 1:n_labels
    edge_endpoints_label_i_filter = find(edge_endpoints_labels == i);
    [ ~, min_ind ] = min(edge_endpoints_pca_filtered(edge_endpoints_label_i_filter, 2));
    [ ~, max_ind ] = max(edge_endpoints_pca_filtered(edge_endpoints_label_i_filter, 2));
    edge_endpoints_filter(edge_endpoints_label_i_filter(min_ind)) = true;
    edge_endpoints_filter(edge_endpoints_label_i_filter(max_ind)) = true;
end
edge_endpoints_filtered = edge_endpoints_filtered(edge_endpoints_filter, :);

if verbose
    figure;
    bw_comparison = cat(3, bw_double, bw_new, zeros(image_height, image_width));
    imshow(bw_comparison);
    hold on
    scatter(edge_endpoints(:, 1), edge_endpoints(:, 2), 'b')
    scatter(edge_endpoints_filtered(:, 1), edge_endpoints_filtered(:, 2), 'g')
    hold off
    title(sprintf(['Original edge image (red) intersected with thresholded filtered edge image (yellow)\n',...
        'Edge endpoints (blue circles) and filtered edge endpoints (green circles) are shown.']))
end

endpoints = edge_endpoints_filtered;

end

