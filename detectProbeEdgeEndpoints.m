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
%   edges in `bw` to join broken edges, and to eliminate portions of edges
%   that protrude along the outside the probe. It is also used to eliminate
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
% and may be broken, due to strong highlights and reflections, but are
% assumed to contain the true points where the edges meet the outer borders
% of the probe.
%
% This function cannot compensate for edges that have been detected shorter
% than they actually are. In general, there is no way to correct for this
% without using information about the physical widths of the junctions
% between bands, which can only be used after the detected edges have been
% matched to probe band junctions.)
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
% - The filtered image is binarized with `threshold` to obtain
%   `bw_filtered`, and intersected with the original edge image to obtain
%   `bw_new`. `bw_new` should contain fewer portions of edges that protrude
%   along the outer borders of the probe.
%
% - Candidate edge endpoints are extracted by finding the extreme points in
%   each binary region in `bw_filtered` along the PCA minor axis of the
%   probe, where the extreme points are constrained to be from the set of
%   white pixels in `bw_new`. The PCA axes used in this step are those
%   computed fresh from the white pixels in `bw_new`.
%
% - Candidate edge endpoints are removed if they are not on the expected
%   side of the PCA major axis. Specifically, the endpoint with the largest
%   second coordinate in the space of the PCA axes should be above the
%   major axis, whereas the other endpoint should be below the major axis.
%   This step will incorrectly reject endpoints if the PCA major axis is
%   located outside the true borders of the probe.
%
% - As there may still be spurious edge fragments lying along the outer
%   borders of the probe, edge endpoints are rejected if their separation
%   is less than `edge_width`. Note that this may incorrectly reject small
%   portions of edges broken by highlights or reflections on the probe.
%   This is not critical as long as these edge portions are not located
%   where the true edges intersect the borders of the probe.
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
        
        % Define the size of the filter so that it covers 1 - 2p of the
        % area under the Gaussian function of distance
        p = 0.05;
        radial_limit = norminv([p (1-p)], 0, s_distance);
        radial_limit = 2 * ceil(abs(radial_limit(2)));
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

% Find edges and their endpoints along the second PCA axis.
% Select candidate points in `bw_new` that have extreme values along the
% second PCA axis. "Extreme" is computed over each edge in `bw_filtered`.
[edge_points_y, edge_points_x] = find(bw_new);
[coeff, score, ~, ~, ~, mu] = pca([edge_points_x, edge_points_y]);
[ bw_filtered_labelled, n_labels ] = bwlabel(bw_filtered);
edge_points_linear = sub2ind(image_size, edge_points_y, edge_points_x);
edge_endpoints_labels = bw_filtered_labelled(edge_points_linear);
edge_endpoints_pca = zeros(n_labels * 2, 1);
edge_endpoints = zeros(n_labels * 2, 2);
k = 1;
for i = 1:n_labels
    edge_endpoints_label_i_filter = (edge_endpoints_labels == i);
    if sum(edge_endpoints_label_i_filter) > 1
        [ edge_endpoints_pca(k), min_ind ] = min(score(edge_endpoints_label_i_filter, 2));
        [ edge_endpoints_pca(k + 1), max_ind ] = max(score(edge_endpoints_label_i_filter, 2));
        edge_points_i = edge_points_linear(edge_endpoints_label_i_filter);
        [ edge_endpoints(k, 2), edge_endpoints(k, 1) ] =...
            ind2sub(image_size, edge_points_i(min_ind));
        [ edge_endpoints(k + 1, 2), edge_endpoints(k + 1, 1) ] =...
            ind2sub(image_size, edge_points_i(max_ind));
        k = k + 2;
    end
end

% Remove endpoints on the same edges that are not separated by the PCA
% major axis of the probe
% Note: `<` and `>` filter out zeros in 'edge_endpoints_pca' corresponding
% to edges for which the test `sum(edge_endpoints_label_i_filter) > 1
% failed above.
edge_endpoints_filter = (edge_endpoints_pca(1:2:end) < 0) &...
    (edge_endpoints_pca(2:2:end) > 0);
edge_endpoints_filter = repelem(edge_endpoints_filter, 2);
edge_endpoints_filtered = edge_endpoints(edge_endpoints_filter, :);

% Remove endpoints on the same edges having a separation on the order of
% the edge width
edge_endpoints_separation = diff(edge_endpoints_filtered, 1, 1);
edge_endpoints_separation = edge_endpoints_separation(1:2:end, :);
edge_endpoints_separation = dot(edge_endpoints_separation, edge_endpoints_separation, 2);
edge_endpoints_filter = edge_endpoints_separation > (edge_width  ^ 2);
edge_endpoints_filter = repelem(edge_endpoints_filter, 2);
edge_endpoints_filtered = edge_endpoints_filtered(edge_endpoints_filter, :);

if verbose
    figure;
    bw_comparison = cat(3, bw_double, bw_new, zeros(image_height, image_width));
    imshow(bw_comparison);
    hold on
    pca_axes = pcaAxes2D( coeff, mu );
    line_points = lineToBorderPoints(pca_axes, image_size);
    line(line_points(1, [1,3])', line_points(1, [2,4])', 'Color', 'c');
    line(line_points(2, [1,3])', line_points(2, [2,4])', 'Color', 'm');
    scatter(edge_endpoints(:, 1), edge_endpoints(:, 2), 'b', 'filled')
    scatter(edge_endpoints_filtered(:, 1), edge_endpoints_filtered(:, 2), 'g', 'filled')
    hold off
    legend('First PCA axis', 'Second PCA axis', 'Original', 'Filtered')
    title(sprintf(['Original edge image (red) intersected with thresholded filtered edge image (yellow)\n',...
        'Edge endpoints (blue circles) and filtered edge endpoints (green circles) are shown.']))
end

endpoints = edge_endpoints_filtered;

end

