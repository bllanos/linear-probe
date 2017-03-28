function [ probe_mask_initial ] = detectBoundingBoxes(...
    I, probe,...
    probe_color_distribution_resolution, probe_color_distributions,...
    probe_color_distribution_increment, rgb_sigma_polyfit,...
    params, varargin...
    )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%% Parse input arguments

nargoutchk(1,1);
narginchk(7,8);

if ~isempty(varargin)
    verbose = varargin{1};
    display_original_image = verbose.display_original_image;
    display_hue_image = verbose.display_hue_image;
    plot_global_hue_estimator = verbose.plot_global_hue_estimator;
    plot_ratio_estimators = verbose.plot_ratio_estimators;
    display_ratio_distribution_backprojections = verbose.display_ratio_distribution_backprojections;
    verbose_region_extraction = verbose.verbose_region_extraction;
    verbose_region_filtering = verbose.verbose_region_filtering;
    display_region_expansion = verbose.display_region_expansion;
else
    display_original_image = false;
    display_hue_image = false;
    plot_global_hue_estimator = false;
    plot_ratio_estimators = false;
    display_ratio_distribution_backprojections = false;
    verbose_region_extraction = false;
    verbose_region_filtering = false;
    display_region_expansion = false;
end

image_width = size(I, 2);
image_height = size(I, 1);
image_n_channels = size(I, 3);
if image_n_channels ~= 3
    error('A 3-channel RGB image is required.')
end
if display_original_image
    figure;
    imshow(I);
    title('Base image');
end

n_colors = size(probe_color_distributions, 2);

%% Compute variables used later

% Obtain hue values
H = rgb2hue(I);

I_double = im2double(I);
R = I_double(:, :, 1);
G = I_double(:, :, 2);
B = I_double(:, :, 3);

if display_hue_image
    figure
    H_color = ones(image_height, image_width, image_n_channels);
    H_color(:, :, 1) = H;
    H_color = hsv2rgb(H_color);
    imshowpair(H, H_color, 'montage');
    title('Hue channel of original image')
end

% Identify which pairs of colours can be adjacent
probe_color_pairs = [probe.colors(1:(end - 1)), probe.colors(2:end)];
probe_color_pairs = sort(probe_color_pairs, 2);
probe_color_pairs = unique(probe_color_pairs, 'rows');

%% Compute the hue variable kernel density estimator for the image

[...
    I_color_distribution,...
    I_color_distribution_increment...
] = hueVariableKernelDensityEstimator(...
    H, R, G, B, true(image_height, image_width),...
    rgb_sigma_polyfit, probe_color_distribution_resolution...
);

if plot_global_hue_estimator
    legend_names = cell(n_colors + 1, 1);
    legend_names{1} = 'Image';
    for i = 1:n_colors
        legend_names{i + 1} = sprintf('Probe colour %d', i);
    end
    plotHueVariableKernelDensityEstimator(...
        I_color_distribution_increment,...
        [I_color_distribution, probe_color_distributions], legend_names...
    );
    title('Hue variable kernel density estimator for the entire image')
end

%% Transform the image using histogram backprojection

% Ratio distributions with respect to the background
ratio_distributions_bg = zeros(...
        probe_color_distribution_resolution, n_colors...
    );
for i = 1:n_colors
    ratio_distributions_bg(:, i) = ratioDistribution(...
            probe_color_distributions(:, i), I_color_distribution...
        );
end

if plot_ratio_estimators
    legend_names = cell(n_colors, 1);
    for i = 1:n_colors
        legend_names{i} = sprintf('Probe colour %d', i);
    end
    plotHueVariableKernelDensityEstimator(...
        probe_color_distribution_increment, ratio_distributions_bg, legend_names...
    );
    title('Ratio hue variable kernel density estimators for probe colours with respect to the background')
end

% Histogram backprojection
ratio_distributions_backprojected_bg = zeros(image_height, image_width, n_colors);
for i = 1:n_colors
    ratio_distributions_backprojected_bg(:, :, i) = queryDiscretized1DFunction(...
            H, ratio_distributions_bg(:, i), I_color_distribution_increment...
        );
end

if display_ratio_distribution_backprojections
    for i = 1:n_colors
        figure
        imshow(...
                ratio_distributions_backprojected_bg(:, :, i) /...
                max(max(ratio_distributions_backprojected_bg(:, :, i)))...
            );
        title(sprintf('Ratio distribution backprojection for probe colour %d with respect to the background', i))
    end
end

%% Find an initial bounding area for the probe

[ probe_regions_initial, probe_regions_bw_initial] = extractBinaryRegions(...
        ratio_distributions_backprojected_bg,...
        params.noise_threshold,...
        params.erosion_radius,...
        [],...
        verbose_region_extraction...
    );

[
    ~,...
    probe_regions_bw_initial_filtered...
] = detectProbeBinaryRegions(...
        probe_regions_initial,...
        probe_regions_bw_initial,...
        probe_color_pairs,...
        params.radius_adj,...
        params.axis_distance_outlier_threshold,...
        verbose_region_filtering...
    );

% Obtain an upper bound on the probe's area
probe_regions_bw_initial_all = any(probe_regions_bw_initial_filtered, 3);
if display_region_expansion
    figure
    imshow(probe_regions_bw_initial_all);
    title('Initial probe colour regions prior to expansion')
end
disk = strel('disk', params.dilation_radius);
probe_regions_bw_initial_all = imdilate(probe_regions_bw_initial_all, disk);
% Find oriented bounding boxes
initial_region_bounds = regionprops(probe_regions_bw_initial_all, {...
    'Centroid', 'Orientation', 'MajorAxisLength', 'MinorAxisLength'...
    });
n_initial_regions = length(initial_region_bounds);
initial_region_rect_corners = cell(n_initial_regions, 1);
for i = 1:n_initial_regions
    half_length = initial_region_bounds(i).MajorAxisLength * params.region_expansion_factor_length / 2;
    half_width = initial_region_bounds(i).MinorAxisLength * params.region_expansion_factor_width / 2;
    angle = initial_region_bounds(i).Orientation;
    length_vector = half_length * [ cosd(angle), -sind(angle) ];
    width_vector = half_width * [ cosd(angle + 90), -sind(angle + 90) ];
    center = initial_region_bounds(i).Centroid;
    initial_region_rect_corners{i} = [
        center + length_vector + width_vector;
        center - length_vector + width_vector;
        center - length_vector - width_vector;
        center + length_vector - width_vector;
        ];
end

% Obtain a mask for the initial bounding boxes
probe_mask_initial = false(image_height, image_width);
for i = 1:n_initial_regions
    probe_mask_initial = probe_mask_initial | roipoly(...
            I, initial_region_rect_corners{i}(:, 1), initial_region_rect_corners{i}(:, 2)...
        );
end
if display_region_expansion
    figure
    probe_regions_bw_initial_all_display = repmat(probe_regions_bw_initial_all, 1, 1, 3);
    probe_regions_bw_initial_all_display(:, :, 1) =...
        probe_regions_bw_initial_all_display(:, :, 1) |...
        imdilate(bwperim(probe_mask_initial), strel('disk', 2));
    imshow(double(probe_regions_bw_initial_all_display));
    title('Initial probe colour regions after dilation, with outline of initial mask')
end
end

