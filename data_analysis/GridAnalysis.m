%% Chequerboard Grid Pose Estimation Accuracy Evaluation

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created January 26, 2018

%% Parameters
input_wildcard = '';
input_files = listFiles(input_wildcard);

camera_params_filename = '';

n_rows = 2;
n_columns = 5;
n_angles = 8;

true_side_length = 30.0;
true_x = [1 3] * true_side_length;
true_y = [-1 0 2 4 6] * true_side_length;
dx = [1 5 5 4 3 2 1 0];
dy = [0 1 3 3 3 3 3 1];

%% Load input data
x = cell(n_rows, n_columns, n_angles);
y = cell(n_rows, n_columns, n_angles);
z = cell(n_rows, n_columns, n_angles);
nx = cell(n_rows, n_columns, n_angles);
ny = cell(n_rows, n_columns, n_angles);
nz = cell(n_rows, n_columns, n_angles);
n_points = cell(n_rows, n_columns, n_angles);

ind = 1;
for i = 1:n_rows
    for j = 1:n_columns
        for k = 1:n_angles
            fileID = fopen(input_files{ind},'r');
            cloud = textscan(fileID, '%f', 'Delimiter', {','});
            cloud = cloud{1};
            cloud = reshape(cloud, 7, []).';
            %cloud = csvread(input_files{ind}); % Does not work for CSV files written on Windows when this code is run on Linux.
            fclose(fileID);
            n_points{i, j, k} = size(cloud, 1);
            x{i, j, k} = cloud(:, 2);
            y{i, j, k} = cloud(:, 3);
            z{i, j, k} = cloud(:, 4);
            nx{i, j, k} = cloud(:, 5);
            ny{i, j, k} = cloud(:, 6);
            nz{i, j, k} = cloud(:, 7);
            ind = ind + 1;
        end
    end
end

angles_deg = atan2d(dy, dx);

load(camera_params_filename, 'cameraParams');
P = [cameraParams.IntrinsicMatrix.' zeros(3, 1)];

%% Align true and measured points

% Project detected points to 2D
[coeff, all_points_pca, ~, ~, ~, mu] = pca([vertcat(x{:}), vertcat(y{:}), vertcat(z{:})]);
% PCA vectors should be along Z, then X, then Y

x_pca = cell(n_rows, n_columns, n_angles);
y_pca = cell(n_rows, n_columns, n_angles);

ind = 1;
end_ind = 1;
for i = 1:n_rows
    for j = 1:n_columns
        for k = 1:n_angles
            end_ind = ind + n_points{i, j, k} - 1;
            x_pca{i, j, k} = all_points_pca(ind:end_ind, 1);
            y_pca{i, j, k} = all_points_pca(ind:end_ind, 2);
            ind = ind + n_points{i, j, k};
        end
    end
end

% Align to the bottom right corner, and align to the bottom edge
origin = mean([x{2, 1, 1}, y{2, 1, 1}, z{2, 1, 1}], 1);
left = mean([x{1, 1, 1}, y{1, 1, 1}, z{1, 1, 1}], 1);
right_vector = origin - left;
right_vector = right_vector ./ repmat(norm(right_vector), 1, 3);
up_vector = coeff(:, 3);
depth_vector = cross(right_vector, up_vector);
depth_vector = depth_vector ./ repmat(norm(depth_vector), 1, 3);
up_vector = cross(depth_vector, right_vector);
up_vector = up_vector ./ repmat(norm(up_vector), 1, 3);
alignment_basis = [right_vector; up_vector; depth_vector].';

true_x = true_x - repelem(true_x(2), n_rows);
true_y = true_y - repelem(true_y(1), n_columns);
[true_x, true_y] = meshgrid(true_x, true_y);
true_x = true_x(:);
true_y = true_y(:);
true_z = zeros(n_rows * n_columns, 1);
true_points = [true_x, true_z, true_y];

true_points_aligned =  (alignment_basis * true_points.').' + repmat(origin, n_rows * n_columns, 1);

%% Plot true and measured points together

% In original coordinates (3D)
for k = [1 n_angles]
    x_k = reshape(x(:, :, k), [], 1);
    x_k = cell2mat(x_k);
    y_k = reshape(y(:, :, k), [], 1);
    y_k = cell2mat(y_k);
    z_k = reshape(z(:, :, k), [], 1);
    z_k = cell2mat(z_k);

    figure;
    hold on
    scatter3(x_k, y_k, z_k, 3, 'k', 'filled')
    scatter3(true_points_aligned(:, 1), true_points_aligned(:, 2), true_points_aligned(:, 3), 'r', 'X')
    hold off
    xlabel('X [mm]')
    ylabel('Y [mm]')
    zlabel('Z [mm]')
    title(sprintf('Points detected at a pointer angle of %g', angles_deg(k)))
    legend('Detected', 'True')
    %axis equal
end

%% Plot RMSE as a function of angle and depth

% Use only one chequerboard column
row = 1;
rmse = cell(n_columns, 1);
depths = true_points_aligned(((row - 1) * n_columns + 1):(row * n_columns), 3);
angles = cell(n_columns, 1);
for j = 1:n_columns
    ind = 1;
    for k = 1:n_angles
        points = [x{row, j, k}, y{row, j, k}, z{row, j, k}];
        n_points = size(points, 1);
        if n_points > 0
            err = points - repmat(true_points_aligned((row - 1) * n_columns + j, :), n_points, 1);
            err = dot(err, err, 2);
            err = sum(err);
            rmse{j}(ind) = sqrt(err / n_points);
            angles{j}(ind) = angles_deg(k);
            ind = ind + 1;
        end
    end
end

line_specs = {'-ok', '-sk', '--dk', '-.^k', ':xk'};
legend_str = cell(n_columns, 1);

figure
hold on
for j = 1:n_columns
    plot(angles{j}, rmse{j}, line_specs{j});
    legend_str{j} = sprintf('Depth %g mm', depths(j));
end
xlabel('Angle on chequerboard (degrees)')
ylabel('RMSE [mm]')
title('RMSE with respect to depth and angle to chequerboard grid')
legend(legend_str)

%% Plot PCA variances as a function of angle and depth

% Use only one chequerboard column again
pca_var = cell(n_columns, 1);
pca_coeff = cell(n_columns, 1);
for j = 1:n_columns
    ind = 1;
    for k = 1:n_angles
        points = [x{row, j, k}, y{row, j, k}, z{row, j, k}];
        n_points = size(points, 1);
        if n_points > 0
            [coeff, ~,latent] = pca(points);
            pca_var{j}(ind, :) = sqrt(latent.');
            pca_coeff{j}(:, :, ind) = coeff;
            ind = ind + 1;
        end
    end
end

dim_names = {'first', 'second', 'third'};

for dim = 1:3

    figure
    hold on
    for j = 1:n_columns
        plot(angles{j}, pca_var{j}(:, dim), line_specs{j});
        legend_str{j} = sprintf('Depth %g mm', depths(j));
    end
    xlabel('Angle on chequerboard (degrees)')
    ylabel('Standard deviation [mm]')
    title(sprintf('Standard deviation along %s principal component', dim_names{dim}))
    legend(legend_str)

end

%% Plot PCA vectors as a function of angle and depth

% for j = 1:n_columns
%
%     vectors = reshape(permute(pca_coeff{j}, [2, 3, 1]), [], 3);
%     scales = reshape(pca_var{j}.', [], 1);
%     vectors = vectors .* repmat(scales, 1, 3);
%     n_coeff = size(pca_coeff{j}, 3);
%     %colors = jet(n_coeff);
%     %colors = [repelem(colors(:, 1), 3), repelem(colors(:, 2), 3), repelem(colors(:, 3), 3)];
%     colors = eye(3);
%     colors = repmat(colors, n_coeff, 1);
%
%     figure;
%     scatter3(vectors(:, 1), vectors(:, 2), vectors(:, 3), [], colors, 'filled');
%
%     % Colours denote first, second, third components
%     xlabel('X [mm]')
%     ylabel('Y [mm]')
%     zlabel('Z [mm]')
%     title(sprintf('PCA vectors at a depth of %g [mm]', true_points_aligned((row - 1) * n_columns + j, 3)))
%     axis equal
% end

% Angles between PCA vectors and probe axis
for j = 1:n_columns
    vectors = reshape(permute(pca_coeff{j}, [2, 3, 1]), [], 3);
    n_coeff = size(pca_coeff{j}, 3);
    probe_axis = zeros(n_coeff, 3);
    for c = 1:n_coeff
        probe_axis(c, :) = mean([nx{row, j, c}, ny{row, j, c}, nz{row, j, c}], 1);
    end
    probe_axis = [repelem(probe_axis(:, 1), 3), repelem(probe_axis(:, 2), 3), repelem(probe_axis(:, 3), 3)];
    dots = dot(vectors, probe_axis, 2);

    figure;
    hold on
    for dim = 1:3
        y = acosd(dots(dim:3:end));
        plot(angles{j}, y, line_specs{dim});
    end
    hold off

    xlabel('Angle on chequerboard (degrees)')
    ylabel('Angle with probe axis (degrees)')
    title(sprintf('PCA vectors at a depth of %g [mm]', true_points_aligned((row - 1) * n_columns + j, 3)))
    legend('First PCA direction', 'Second PCA direction', 'Third PCA direction')
end

% Look at PCA vectors in the image plane
for j = 1:n_columns
    vectors = reshape(permute(pca_coeff{j}, [2, 3, 1]), [], 3);
    n_coeff = size(pca_coeff{j}, 3);
    probe_axis = zeros(n_coeff, 3);
    for c = 1:n_coeff
        probe_axis(c, :) = mean([nx{row, j, c}, ny{row, j, c}, nz{row, j, c}], 1);
    end

    % Project to image space
    vectors = [vectors, zeros(size(vectors, 1), 1)]; %#ok<AGROW>
    probe_axis = [probe_axis, zeros(size(probe_axis, 1), 1)]; %#ok<AGROW>
    image_vectors = (P * vectors.').';
    image_vectors = image_vectors ./ repmat(image_vectors(:, end), 1, 3);
    image_vectors = image_vectors(:, 1:2);
    image_vectors = image_vectors ./ sqrt(dot(image_vectors, image_vectors, 2));
    image_probe_axis = (P * probe_axis.').';
    image_probe_axis = image_probe_axis ./ repmat(image_probe_axis(:, end), 1, 3);
    image_probe_axis = [repelem(image_probe_axis(:, 1), 3), repelem(image_probe_axis(:, 2), 3)];
    image_probe_axis = image_probe_axis ./ sqrt(dot(image_probe_axis, image_probe_axis, 2));
    dots = dot(image_vectors, image_probe_axis, 2);

    figure;
    hold on
    for dim = 2:3
        y = acosd(dots(dim:3:end));
        plot(angles{j}, y, line_specs{dim});
    end
    hold off

    xlabel('Angle on chequerboard (degrees)')
    ylabel('Angle with probe axis in image (degrees)')
    title(sprintf('PCA vectors at a depth of %g [mm]', true_points_aligned((row - 1) * n_columns + j, 3)))
    legend('Second PCA direction', 'Third PCA direction')
end
