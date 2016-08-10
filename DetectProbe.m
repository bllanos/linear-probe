%% Probe Object Detection
% Use previously estimated colour distributions and measured dimensions of
% the probe object to locate it in a test image.
%
% ## Usage
%   Modify parameters and paths to input data in the first code section
%   below, then run.
%
% ## Probe Detection Assumptions and Limitations
% - The image only contains a single instance of the probe.
% - The probe may be occluded in one or more places. However, it can be
%   detected if at least 5 transition lines between pairs of coloured bands
%   of the probe are completely visible. The transition lines do not need
%   to be between bands that are all consecutive.
%
% ## Input
%
% ### Probe detection model
% A '.mat' file containing several variables, which is the output of
% '.\CreateProbeDetectionModel.m'. Refer to the documentation of
% '.\CreateProbeDetectionModel.m' for details.
%
% ### Image containing the probe
% An RGB image in any format that can be loaded by `imread` showing enough
% of the probe for detection to be possible. The image should have been
% taken by the same camera, and under the same camera parameters, as the
% image used to create the probe detection model. Ideally, the lighting
% conditions in the image are similar to those in the context of which the
% detection model was created.
%
% ### Colour noise parameters
% A '.mat' file containing a 'rgb_sigma_polyfit' variable, as output by the
% script '.\EstimateRGBStandardDeviations.m'. 'rgb_sigma_polyfit' describes
% the variation in RGB channel standard deviations with RGB values in the
% image. This information should be computed from images taken
% under the same conditions and with the same camera parameters as the
% image in which the probe is to be detected, if not computed from this
% image.
%
% ## Output
%
% ## References
% - M.-C. Chuang, J.-N. Hwang, K. Williams and R. Towler. "Tracking Live
%   Fish from Low-Contrast and Low-Frame-Rate Stereo Videos". IEEE
%   Transactions on Circuits and Systems for Video Technology, vol. 25, no.
%   1, pp. 167-179, Jan. 2015.
% - T. Gevers and H. Stokman. "Robust Histogram Construction from Color
%   Invariants for Object Recognition". IEEE Transactions on Pattern
%   Analysis and Machine Intelligence, vol. 26, no. 1, pp. 113-118, Jan.
%   2004.

% Bernard Llanos
% Spring 2016 research assistantship supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created August 10, 2016

%% Input data and parameters

% List of parameters to save with results
parameters_list = {
    };

% Probe detection model
detection_model_filename = 'C:\Users\Bernard\Documents\Data\20160725_probes_lotus\20160810_redPen_probeDetectionModel_bottomCamera.mat';
% Image of probe in an unknown pose
I_filename = 'C:\Users\Bernard\Documents\Data\20160725_probes_lotus\probesAgainstWhiteBox\original\redPen_white_b_1.bmp';
% RGB noise parameters
rgb_sigma_filename = 'C:\Users\Bernard\Documents\Data\20160725_probes_lotus\20160809_rgbStddev_redPen_bottomCamera.mat';

% Debugging tools
display_original_image = false;
display_hue_image = false;
plot_global_hue_estimator = false;
plot_ratio_estimators = false;
display_ratio_distribution_backprojections = true;

%% Load the image containing the probe in an unknown pose

I = imread(I_filename);
image_width = size(I, 2);
image_height = size(I, 1);
image_n_channels = size(I, 3);
if image_n_channels ~= 3
    error('A 3-channel RGB image is required.')
end
if display_original_image
    figure; %#ok<UNRCH>
    imshow(I);
    title('Base image');
end

%% Load the probe detection model
model_variables_required = {...
        'probe_band_color_distribution_resolution',...
        'probe',...
        'probe_band_color_distributions',...
        'probe_band_color_distribution_increment'...
    };
load(detection_model_filename, model_variables_required{:});
if ~all(ismember(model_variables_required, who))
    error('One or more of the probe detection model variables is not loaded.')
end

%% Compute the hue variable kernel density estimator for the image

% Obtain hue values
H = rgb2hue(I);
    
if display_hue_image
    figure %#ok<UNRCH>
    imshow(H);
    title('Hue channel of original image')
end

% Compute hue variable kernel density estimators from probe bands
load(rgb_sigma_filename, 'rgb_sigma_polyfit');
if ~exist('rgb_sigma_polyfit', 'var')
    error('No variable called ''rgb_sigma_polyfit'' is loaded (which would contain the camera RGB noise model).')
end

I_double = im2double(I);
R = I_double(:, :, 1);
G = I_double(:, :, 2);
B = I_double(:, :, 3);
[...
    I_color_distribution,...
    I_color_distribution_increment...
] = hueVariableKernelDensityEstimator(...
    H, R, G, B, true(image_height, image_width),...
    rgb_sigma_polyfit, probe_band_color_distribution_resolution...
);

if plot_global_hue_estimator
    x = 0:I_color_distribution_increment:1; %#ok<UNRCH>
    figure
    plot(x, I_color_distribution)
    title('Hue variable kernel density estimators for the entire image')
    xlabel('Hue, \theta (range [0, 1])')
    ylabel('Density, P(\theta)')
end

%% Transform the image using histogram backprojection
n_bands = size(probe_band_color_distributions, 2);

ratio_color_distributions = zeros(...
        probe_band_color_distribution_resolution, n_bands...
    );
for i = 1:n_bands
    ratio_color_distributions(:, i) = ratioDistribution(...
            probe_band_color_distributions(:, i), I_color_distribution...
        );
end

if plot_ratio_estimators
    x = 0:probe_band_color_distribution_increment:1; %#ok<UNRCH>
    line_styles = {'-', '--', ':', '-.'};
    legend_names = cell(n_bands, 1);
    plot_colors = jet(n_bands);
    figure
    hold on
    for i = 1:n_bands
        legend_names{i} = sprintf('Probe band %d', i);
        plot(...
                x, ratio_color_distributions(:, i),...
                'Color', plot_colors(i, :),...
                'LineStyle', line_styles{mod(i - 1, length(line_styles)) + 1}...
            )
    end
    hold off
    legend(legend_names{:});
    title('Ratio hue variable kernel density estimators for bands on the probe')
    xlabel('Hue, \theta (range [0, 1])')
    ylabel('Density, P(\theta)')
end

% Histogram backprojection
ratio_color_distributions_backprojected = zeros(image_height, image_width, n_bands);
for i = 1:n_bands
    ratio_color_distributions_backprojected(:, :, i) = queryDiscretized1DFunction(...
            H, ratio_color_distributions(:, i), I_color_distribution_increment...
        );
end

if display_ratio_distribution_backprojections
    for i = 1:n_bands %#ok<UNRCH>
        figure
        imshow(...
                ratio_color_distributions_backprojected(:, :, i) /...
                max(max(ratio_color_distributions_backprojected(:, :, i)))...
            );
        title(sprintf('Ratio distribution backprojection for probe band %d', i))
    end
end