%% RGB Noise Estimation (Non-Interactive)
% Create a rough model of RGB standard deviations from several images or a
% video. The images or video frames are assumed to show the same scene, and
% differ only because of image noise.
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
% ### RGB noise model
% A '.mat' file containing the following variables:
% - 'rgb_sigma_polyfit': A 2 x n_channels array, containing the coefficients
%   of lines fitted to `(mean(channel_values),var(channel_values))` pairs.
%   'n_channels' is 3 (for Red, Green, Blue, in that order).
%   'rgb_sigma_polyfit' describes the relationship between colour channel
%   values and colour channel noise, as calculated using MATLAB's `polyfit`
%   function. Specifically, each column of 'rgb_sigma_polyfit' corresponds
%   to the first output argument, 'p', of `polyfit`.
%
%   Note that colour channel values are in the range [0, 1] (and standard
%   deviations correspond to this range as well, therefore).
%
% - If a video file was used as input:
%   - 'video_filename': A string containing the filepath of the input video
%   - 'start_time': The time offset in the video at which to start
%     collecting frames for estimating image noise.
%   - 'n_frames': The number of frames in the video used to estimate image
%     noise.
%
% - If a directory of image files was used as input:
%   - 'in_directory': A cell vector containing the directory holding the
%     input images
%   - 'wildcard': The wildcard used to select image files using `ls()`
%
% ## References
% - Code created in August 2017 for loading RAW images.
% - 'EstimateRGBStandardDeviationsInteractive.m'
% - Martinec, E. (2008). Noise, dynamic range and bit depth in digital SLR.
%   Retrieved from http://theory.uchicago.edu/∼ejm/pix/20d/tests/noise/

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created January 8, 2018

%% Input data and parameters

% Process a batch of images, or a video
video_mode = true;

if video_mode
    % Input video file
    video_filename = '/home/llanos/GoogleDrive/PointProbing/DataAndResults/20180212_reverseEXPO/noise/noise_csc265.webm';
    
    % Starting offset in the video
    start_time = 0;
    
    % Maximum number of frames to read (set to zero to read all frames)
    n_frames = 0;
else
    % Directory containing the input images
    in_directory = '/home/llanos/GoogleDrive/PointProbing/DataAndResults/20180108_noiseEstimation';
    
    % Input image filename wildcard
    wildcard = '*.jpg';
end

% Directory in which to save the output file
out_directory = '/home/llanos/GoogleDrive/PointProbing/DataAndResults/20180212_reverseEXPO/noise';

% Name of the output file
out_filename = 'noise_csc265_rgbstddev_nonInteractive.mat';

% Number of pixels to plot when visualizing the results
plot_count = 1000;

%% Process the images
if video_mode
    frames = readVideo(video_filename, start_time, n_frames);
    n_frames = size(frames, 4);
else
    frames = readImageGroup(in_directory, wildcard);
end
n_channels = size(frames, 3);

px_mean = mean(frames, 4);
px_var = var(frames, 0, 4);
px_mean = reshape(px_mean, [], n_channels);
px_var = reshape(px_var, [], n_channels);
n_px = size(px_mean, 1);

% Fit a line to the square of the image noise vs. the pixel intensities
rgb_sigma_polyfit = zeros(2, n_channels);
for i = 1:n_channels
    rgb_sigma_polyfit(:, i) = polyfit(px_mean(:, i), px_var(:, i), 1);
end

% Plot the results
figure;
hold on
markerspec = {'ro', 'go', 'bo'};
linespec = {'r--', 'g--', 'b--'};
x = linspace(0, 1, 50);
plot_sample = [px_mean, px_var];
plot_sample = plot_sample(randperm(n_px,plot_count), :);
for i = 1:n_channels
    scatter(plot_sample(:, i), sqrt(plot_sample(:, i + n_channels)), markerspec{i});
    plot(x, sqrt(polyval(rgb_sigma_polyfit(:, i), x)), linespec{i});
end
hold off
title('RGB standard deviations')
xlabel('Mean colour value')
ylabel('Standard deviation')

%% Save the preprocessed data
if video_mode
    save(...
        fullfile(out_directory, out_filename),...
        'rgb_sigma_polyfit', 'video_filename', 'start_time', 'n_frames'...
    );
else
    save(...
        fullfile(out_directory, out_filename),...
        'rgb_sigma_polyfit', 'in_directory', 'wildcard'...
    );
end