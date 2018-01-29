%% 3D reconstruction of the square bottom of a container for paper clips

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created January 28, 2018

%% Data input parameters
input_file = '/home/llanos/GoogleDrive/PointProbing/DataAndResults/20180126_paperClipCase/case_square/tour.csv';

% start_analysis = 621;
% end_analysis = 719;
start_analysis = 1;
end_analysis = 719;

%% Load input data

fileID = fopen(input_file,'r');
cloud = textscan(fileID, '%f', 'Delimiter', {','});
cloud = cloud{1};
cloud = reshape(cloud, 7, []).';
%cloud = csvread(input_files{ind}); % Does not work for CSV files written on Windows when this code is run on Linux.
fclose(fileID);
n_points_all = size(cloud, 1);
x_all = cloud(:, 2);
y_all = cloud(:, 3);
z_all = cloud(:, 4);
nx_all = cloud(:, 5);
ny_all = cloud(:, 6);
nz_all = cloud(:, 7);

x = x_all(start_analysis:end_analysis);
y = y_all(start_analysis:end_analysis);
z = z_all(start_analysis:end_analysis);
nx = nx_all(start_analysis:end_analysis);
ny = ny_all(start_analysis:end_analysis);
nz = nz_all(start_analysis:end_analysis);
n_points = end_analysis - start_analysis + 1;

%% Visualize the points

point_colors_all = 1:n_points_all;

% % All points
% figure;
% hold on
% scatter3(x_all, y_all, z_all, [], point_colors_all, 'filled');
% %quiver3(x, y, z, nx, nz, ny, 'Color', 'black');
% hold off
% colormap jet
% c = colorbar;
% c.Label.String = 'Frame number';
% 
% xlabel('X')
% ylabel('Y')
% zlabel('Z')
% axis equal
% grid on
% title('Points on tours of the square')
% 
% % Points of interest
% point_colors = 1:n_points;
% figure;
% hold on
% scatter3(x, y, z, [], point_colors, 'filled');
% %quiver3(x, y, z, nx, nz, ny, 'Color', 'black');
% hold off
% colormap jet
% c = colorbar;
% c.Label.String = 'Frame number';
% 
% xlabel('X')
% ylabel('Y')
% zlabel('Z')
% axis equal
% grid on
% title('Points from the portion of interest')

%% Select a bounding box

% Bounding box
x_lim = [-147 -93];
y_lim = [-104 -50];
z_lim = [410 502];

cloud_box = [x y z nx ny nz];
cloud_box = cloud_box(...
    cloud_box(:, 1) > x_lim(1) & cloud_box(:, 1) < x_lim(2) &...
    cloud_box(:, 2) > y_lim(1) & cloud_box(:, 2) < y_lim(2) &...
    cloud_box(:, 3) > z_lim(1) & cloud_box(:, 3) < z_lim(2)...
    , :);

% Points of interest within the box
figure;
hold on
scatter3(cloud_box(:, 1), cloud_box(:, 2), cloud_box(:, 3), 'filled');
%quiver3(x, y, z, nx, nz, ny, 'Color', 'black');
hold off

xlabel('X')
ylabel('Y')
zlabel('Z')
axis equal
grid on
title('Points from the bounding box within the portion of interest')

%% Express in PCA space

[coeff, ~, ~, ~, ~, mu] = pca(cloud_box(:, 1:3));
transform_pca = [coeff, mu.'; 0 0 0 1];
cloud_box_pca =  (transform_pca \ [cloud_box(:, 1:3), ones(size(cloud_box, 1), 1)].').';

% figure;
% scatter3(cloud_box_pca(:, 1), cloud_box_pca(:, 2), cloud_box_pca(:, 3), 'filled');
% xlabel('First PCA component')
% ylabel('Second PCA component')
% zlabel('Third PCA component')
% axis equal
% grid on
% title('Points from the bounding box within the portion of interest (PCA space)')

% Points manually selected from the plot
corners_pca = [
    -4.084 -22.56 0.9063;
    26.58 -1.338 4.195;
    -1.036 27.35 1.803;
    -27.11 5.974 -4.331
    ];
n_corners = size(corners_pca, 1);
corners = (transform_pca * [corners_pca ones(n_corners, 1)].').';

% Filter to points within a distance from lines joining the corners
distance_threshold = 2;
distances = zeros(size(cloud_box_pca, 1), n_corners);
for c = 1:n_corners
    endpoints = corners_pca([c (mod(c, n_corners) + 1)], :);
    distances(:, c) = distanceToLine(cloud_box_pca(:, 1:3), endpoints);
end
distance_filter = any(distances <= distance_threshold, 2);

square_pca = cloud_box_pca(distance_filter, :);

figure;
scatter3(square_pca(:, 1), square_pca(:, 2), square_pca(:, 3), 'filled');
xlabel('First PCA component')
ylabel('Second PCA component')
zlabel('Third PCA component')
axis equal
grid on
title('Filtered points for ICP')

%% Project onto an approximate plane of the square

square = (transform_pca * square_pca.').';
[coeff, ~, ~, ~, ~, mu] = pca(square(:, 1:3));
transform_pca = [coeff, mu.'; 0 0 0 1];
cloud_box_pca =  (transform_pca \ [cloud_box(:, 1:3), ones(size(cloud_box, 1), 1)].').';

figure;
scatter3(cloud_box_pca(:, 1), cloud_box_pca(:, 2), cloud_box_pca(:, 3), 'filled');
xlabel('First PCA component')
ylabel('Second PCA component')
zlabel('Third PCA component')
axis equal
grid on
title('Points from the bounding box within the portion of interest (Filtered PCA space)')

%% Register with the ideal square using ICP

square_side_length = 36.5;
edge_sampling = 20;
square_true = [
    linspace(0, square_side_length, edge_sampling).' zeros(edge_sampling, 1);
    repmat(square_side_length, edge_sampling, 1) linspace(0, square_side_length, edge_sampling).';
    linspace(square_side_length, 0, edge_sampling).' repmat(square_side_length, edge_sampling, 1);
    zeros(edge_sampling, 1) linspace(square_side_length, 0, edge_sampling).'
    ];

% Give an initial approximate alignment
corners_basis_1 = corners(2, 1:3) - corners(1, 1:3);
corners_basis_1 = corners_basis_1 ./ repmat(norm(corners_basis_1), 1, 3);
corners_basis_2 = corners(4, 1:3) - corners(1, 1:3);
corners_basis_2 = corners_basis_2 ./ repmat(norm(corners_basis_2), 1, 3);
corners_basis_3 = cross(corners_basis_1, corners_basis_2);
corners_basis_3 = corners_basis_3 ./ repmat(norm(corners_basis_3), 1, 3);
corners_basis_2 = cross(corners_basis_3, corners_basis_1);
corners_transform = [corners_basis_1.' corners_basis_2.' corners_basis_3.' corners(1, 1:3).'; 0 0 0 1];

square_true_data_space = (corners_transform * [square_true zeros(edge_sampling * n_corners, 1) ones(edge_sampling * n_corners, 1)].').';

% ICP
fixed = pointCloud(square(:, 1:3));
moving = pointCloud(square_true_data_space(:, 1:3));

icp_tform = pcregrigid(moving,fixed);
icp_tform = icp_tform.T';
registration_tform = icp_tform * corners_transform;

% Display the result in the original data space
square_true_corners = [
    0 0;
    square_side_length 0;
    square_side_length square_side_length;
    0 square_side_length
];
square_true_corners_data_space = (registration_tform * [square_true_corners zeros(n_corners, 1) ones(n_corners, 1)].').';

figure;
scatter3(x, y, z, 'k', '.');
hold on
for c = 1:n_corners
    endpoints = square_true_corners_data_space([c (mod(c, n_corners) + 1)], 1:3);
    line(endpoints(:, 1),endpoints(:, 2),endpoints(:, 3),'Color','red','LineWidth',2)
end
hold off
xlabel('X')
ylabel('Y')
zlabel('Z')
axis equal
grid on
title('Datapoints and registered square (in camera space)')

% Display the result in the space of the true corners
points_canonical_space = (registration_tform \ [x y z ones(n_points, 1)].').';

figure;
scatter3(points_canonical_space(:, 1), points_canonical_space(:, 2), points_canonical_space(:, 3), 'k', '.');
hold on
for c = 1:n_corners
    endpoints = [
        square_true_corners([c (mod(c, n_corners) + 1)], :) zeros(2, 1)
        ];
    line(endpoints(:, 1),endpoints(:, 2),endpoints(:, 3),'Color','red','LineWidth',2)
end
hold off
xlabel('X')
ylabel('Y')
zlabel('Z')
axis equal
grid on
title('Datapoints and registered square (in canonical space)')