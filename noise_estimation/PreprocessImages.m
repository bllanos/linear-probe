%% Preprocess Images
% - Read and average each group of images, and load their pixel class labels.
% - Fit a model for image noise
%
% ## Usage
% Modify the parameters, the first code section below, then run.
%
% ## Input
%
% Refer to the first code section below.
%
% ## Output
%
% A '.mat' file containing the preprocessed image data.
%
% ## References
% - Code created in August 2017 for loading RAW images.
% - Code created in August 2016 for modelling image noise
% - Martinec, E. (2008). Noise, dynamic range and bit depth in digital SLR.
%   Retrieved from http://theory.uchicago.edu/∼ejm/pix/20d/tests/noise/

% Created for: CMPUT 551 Mini project: Colour classification
% Fall 2017
% Bernard Llanos, ID 1505236
% Department of Computing Science, University of Alberta

%% Input data and parameters

% Directory containing the input images
in_directory = '/home/llanos/Downloads/data';

% Directory in which to save the output images
out_directory = './dataset';

% Input filename wildcards
% Filenames must not contain spaces
wildcards = {
    'tape_fluorescentLight_2*',...
    'tape_noLight_2*',...
    'tape_workLight_2*',...
    'tape_workLight_fluorescentLight_2*'
};

% Number of pixels to plot when visualizing the dataset
plot_count = 1000;

%% Process the images
n_groups = length(wildcards);
[px_1, px_std_1, sz] = loadImageGroup(in_directory, wildcards{1});
n_px = size(px_1, 1);
n_channels = size(px_1, 2);
px = zeros(n_px, n_channels, n_groups);
px(:, :, 1) = px_1;
px_std = zeros(n_px, size(px_std_1, 2), n_groups);
px_std(:, :, 1) = px_std_1;
for i = 2:n_groups
    [px(:, :, i), px_std(:, :, i)] = loadImageGroup(in_directory, wildcards{i});
end

%% Load class labels
I_labels = imread(fullfile(in_directory, 'labels.bmp'));
I_labels = I_labels(:, :, 1);
I_labels = reshape(I_labels, [], 1);
n_classes = max(I_labels);
class_sizes = zeros(n_classes, 1); % Class sizes are the same across lighting conditions
class_filters = false(n_px, n_classes);
for i = 1:n_classes
    class_filters(:, i) = (I_labels == i);
    class_sizes(i) = sum(class_filters(:, i));
end

% Group pixels into classes: Pixels organized by class, then by lighting
% condition.
px_labelled_classFirst = cell(n_classes, 1);
for i = 1:n_classes
    px_labelled_classFirst{i} = cell(n_groups, 1);
    for j = 1:n_groups
        px_labelled_classFirst{i}{j} = px(class_filters(:, i), :, j);
    end
end

%% Produce alternate organizations of the data

% Pixels organized by lighting condition, then by class
px_labelled_groupFirst = cell(n_groups, 1);
for j = 1:n_groups
    px_labelled_groupFirst{j} = cell(n_classes, 1);
    for i = 1:n_classes
        px_labelled_groupFirst{j}{i} = px_labelled_classFirst{i}{j};
    end
end

% Pixels organized by class, with lighting conditions combined
px_labelled_classOnly = cell(n_classes, 1);
for i = 1:n_classes
    px_labelled_classOnly{i} = vertcat(px_labelled_classFirst{i}{:});
end

%% Visualize the dataset

class_names = {
    'black', 'white', 'pink', 'red', 'orange', 'yellow', 'green', 'sky blue', 'deep blue'...
    }.';

lighting_names = {
    'fluorescent light',...
    'no light',...
    'work light',...
    'both lights'
};

% Pixels organized by class, then by lighting condition
% One plot per class
plotRGBGroups(px_labelled_classFirst, plot_count, lighting_names, class_names);

% Pixels organized by lighting condition, then by class
% One plot per lighting condition
plotRGBGroups(px_labelled_groupFirst, plot_count, class_names, lighting_names);

% Pixels organized by class, with lighting conditions combined
% One plot
plotRGB(px_labelled_classOnly, plot_count, class_names)
title('Colour classes over all lighting conditions')

%% Model image noise

% Fit a line to the square of the image noise vs. the pixel intensities
polyfit_degree = 1;
rgb_sigma_polyfit = zeros(polyfit_degree + 1, n_channels);
px_std_all = reshape(permute(px_std, [3 1 2]), [], n_channels);
px_std_all_sq = px_std_all .^ 2;
px_all = reshape(permute(px, [3 1 2]), [], n_channels);
for i = 1:n_channels
    rgb_sigma_polyfit(:, i) = polyfit(px_all(:, i), px_std_all_sq(:, i), polyfit_degree);
end

% Plot the results
figure;
hold on
markerspec = {'ro', 'go', 'bo'};
linespec = {'r--', 'g--', 'b--'};
x = linspace(0, 1, 50);
n_samples = 1e2;
plot_sample = datasample([px_all, px_std_all], n_samples);
for i = 1:n_channels
    scatter(plot_sample(:, i), plot_sample(:, i + n_channels), markerspec{i});
    plot(x, sqrt(polyval(rgb_sigma_polyfit(:, i), x)), linespec{i});
end
hold off
title('RGB standard deviations')
xlabel('Mean colour value')
ylabel('Standard deviation')

%% Save the preprocessed data
save(...
    fullfile(out_directory, 'data.mat'),...
    'class_names', 'class_sizes', 'lighting_names', 'n_channels',...
    'n_classes', 'n_groups', 'polyfit_degree', 'px_labelled_classFirst',...
    'px_labelled_classOnly', 'px_labelled_groupFirst', 'rgb_sigma_polyfit',...
    'sz'...
    )