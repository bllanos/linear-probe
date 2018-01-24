%% Probe Object Detection
% Use previously estimated colour distributions and measured dimensions of
% the probe object to locate it in a test image.
%
% ## Usage
%   Modify the parameters and the paths to input data in the first code
%   section below, then run.
%
% ## Probe Detection Assumptions and Limitations
% - The image contains either saturated colours which are different from
%   the colours of the probe, or fairly unsaturated colours, to assist with
%   colour discrimination using hues and thresholds on colour saturation.
% - The image only contains a single instance of the probe.
% - The probe may be occluded in one or more places. However, it can be
%   located in 3D space if at least 3 edges between pairs of coloured bands
%   of the probe are completely visible. (Localization is not handled by
%   this script, but by 'LocalizeProbeScript.m'.) The edges do not need to
%   be between bands that are all consecutive.
%
% ## Input
%
% ### Probe detection model
% A '.mat' file containing several variables, which is the output of
% 'CreateProbeDetectionModel.m'. Refer to the documentation of
% 'CreateProbeDetectionModel.m' for details.
%
% ### Image containing the probe
% An RGB image, in any format that can be loaded by `imread`, showing enough
% of the probe for detection to be possible. The image has ideally been
% taken by the same camera, and under the same camera parameters, as the
% image used to create the probe detection model. Ideally, the lighting
% conditions in the image are similar to those in the context of which the
% detection model was created.
%
% The image should have been corrected for lens distortion.
%
% ### Colour noise parameters
% Required if at least one of the parameters 'uniform_background_initial',
% or 'uniform_background_final' is false.
%
% A '.mat' file containing a 'rgb_sigma_polyfit' variable, as output by the
% script '.\EstimateRGBStandardDeviations.m'. 'rgb_sigma_polyfit' describes
% the variation in RGB channel standard deviations with RGB values in the
% image. This information should be computed from images taken
% under the same conditions and with the same camera parameters as the
% image in which the probe is to be detected, if not computed from this
% same image.
%
% 'rgb_sigma_polyfit' is used to compute a colour distribution for the
% image background. Instead, if 'uniform_background_initial' is true, a
% uniform distribution is used during the first pass. Likewise, if
% 'uniform_background_final' is true, a uniform distribution is used during
% the second pass.
%
% ## Output
%
% ### Probe detection results
% A '.mat' file containing the following variables:
%
% - 'model_filename': The path to the file containing user-provided
%   measurements of the probe in the structure 'probe'. The 'probe'
%   structure is provided to this script via the output of
%   'CreateProbeDetectionModel.m', and so 'model_filename' is actually a
%   direct parameter of 'CreateProbeDetectionModel.m' and an indirect
%   parameter of this script. It is retrieved from the output of
%   'CreateProbeDetectionModel.m', and copied to the output of this script
%   for completeness.
%
% - 'probe_detection_matches': The `matches` output argument of
%   'detectProbe()'. Refer to the documentation of 'detectProbe.m' for
%   details.
%
% - 'probe_detection_matches_filtered': A copy of
%   'probe_detection_matches', with detected edges that were not matched to
%   probe measurements filtered out. The `matches_filtered` output argument
%   of 'detectProbe()'. Refer to the documentation of 'detectProbe.m' for
%   details.
%
% Additionally, the output file contains the values of all parameters in
% the first section of the script below, for reference. (Specifically,
% those listed in `parameters_list`, which should be updated if the set of
% parameters is changed.)
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
        'detection_model_filename',...
        'I_filename',...
        'rgb_sigma_filename',...
        'use_kernel_estimators',...
        'detectionParams'...
    };

% Probe detection model
detection_model_filename = 'C:\Users\llanos\Google Drive\PointProbing\DataAndResults\20180123_bluePenGreenTape_flea3\model\bluePenGreenTape_detectionModel.mat';
% Image of probe in an unknown pose
I_filename = 'C:\Users\llanos\Google Drive\PointProbing\DataAndResults\20180123_bluePenGreenTape_flea3\model\bluePenGreenTape.bmp';
% RGB noise parameters
rgb_sigma_filename = 'C:\Users\llanos\Google Drive\PointProbing\DataAndResults\20180123_bluePenGreenTape_flea3\model\rgbstddev_nonInteractive_video.mat';

% Parameters which do not usually need to be changed
run('SetFixedParameters.m')

%% Load the image containing the probe in an unknown pose
I = imread(I_filename);

%% Load the probe detection model
model_variables_required = { 'probe', 'model_filename' };
if use_kernel_estimators
    model_variables_required(end + 1) = {'probe_color_distributions_kernel'};
else
    model_variables_required(end + 1) = {'probe_color_distributions_gaussian'};
end
load(detection_model_filename, model_variables_required{:});
if ~all(ismember(model_variables_required, who))
    error('One or more of the probe detection model variables is not loaded.')
end
if use_kernel_estimators
    probe_color_distributions = probe_color_distributions_kernel;
else
    probe_color_distributions = probe_color_distributions_gaussian;
end

if ~uniform_background_initial || ~uniform_background_final
    load(rgb_sigma_filename, 'rgb_sigma_polyfit');
    if ~exist('rgb_sigma_polyfit', 'var')
        error('No variable called ''rgb_sigma_polyfit'' is loaded (which would contain the camera RGB noise model).')
    end
else
    rgb_sigma_polyfit = [];
end

%% Detection

[ probe_detection_matches_filtered, probe_detection_matches ] = detectProbe(...
       I, probe, probe_color_distributions, rgb_sigma_polyfit, detectionParams, verbose...
);

%% Save results to a file
save_variables_list = [ parameters_list, {...
        'model_filename',...
        'probe_detection_matches',...
        'probe_detection_matches_filtered'
    } ];
uisave(save_variables_list,'probeDetectionResults');