function [ endpoints ] = detectProbeEdgeEndpoints( bw, edge_width, sigma_theta, threshold, varargin )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%
% varargin = verbose

    function [ H ] = directionalKernel(s_distance, s_theta, mu_theta)
        % Returns a filtering kernel which responds to pixels pixels along
        % lines of high intensity with orientations close to `mu_theta`.
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
edge_endpoints_filter = edge_endpoints_separation > edge_width;
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

