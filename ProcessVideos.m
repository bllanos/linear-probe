%% Track Probe in Video
% Detect and localize the probe in video frames.
%
% ## Usage
%   Modify the parameters and the paths to input data in the first code
%   section below, then run.
%
% ## Probe Detection Assumptions and Limitations
% - See notes on this topic in 'DetectProbeScript.m'.
%
% ## Input
%
% ### Probe detection model
% A '.mat' file containing several variables, which is the output of
% 'CreateProbeDetectionModel.m'. Refer to the documentation of
% 'CreateProbeDetectionModel.m' for details.
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
% videos in which the probe is to be detected.
%
% 'rgb_sigma_polyfit' is used to compute a colour distribution for the
% image background. Instead, if 'uniform_background_initial' is true, a
% uniform distribution is used during the first colour detection pass.
% Likewise, if 'uniform_background_final' is true, a uniform distribution
% is used during the second colour pass.
%
% ### Camera calibration
% A '.mat' file containing a 'cameraParams' variable. The variable is a
% structure of class 'cameraParameters' describing the camera, such as
% output by the MATLAB Camera Calibrator app.
%
% Note: The camera should be calibrated using the same units of measurement
% as used for the measurements of the probe in the file referred to by
% the 'model_filename' variable in the file referred to by
% 'detection_model_filename' below.
%
% ### Videos containing the probe
%
% This script can process either live video from a webcam, or one or more
% saved videos. Refer to the MATLAB documentation for the
% 'vision.VideoFileReader' system object concerning compatible video file
% formats. The videos should have been captured under the same camera
% parameters and, ideally, the same lighting conditions, as used during
% probe detection model creation.
%
% Video frames will be corrected for lens distortion during processing.
%
% ## Output
%
% ### Probe detection and localization results
% '.mat' files containing the following variables:
%
% - 'model_filename': The path to the file containing user-provided
%   measurements of the probe in the structure 'probe'. The 'probe'
%   structure is provided to this script via the output of
%   'CreateProbeDetectionModel.m', and so 'model_filename' is a direct
%   parameter of 'CreateProbeDetectionModel.m' and an indirect parameter of
%   this script. It is retrieved from the output of
%   'CreateProbeDetectionModel.m', and copied to the output of this script
%   for completeness.
%
% - 'localizations': The `localizations` output argument of
%   'trackInVideo()'. Refer to the documentation of 'trackInVideo.m' for
%   details.
%
% - 'detections': The `detections` output argument of 'trackInVideo()'.
%   Refer to the documentation of 'trackInVideo.m' for details.
%
% Additionally, the output files contain the values of all parameters in
% the first section of the script below, for reference. (Specifically,
% those listed in `parameters_list`, which should be updated if the set of
% parameters is changed.)
%
% One output file will be produced per input video, and saved in the
% directory referred to by 'output_data_directory' below. If
% 'output_data_directory' is empty, no output data files will be produced.
%
% ### Annotated videos
%
% The input videos will be annotated with the reprojection of the probe,
% and saved in the directory referred to by 'output_video_directory' below.
% If 'output_video_directory' is empty, no output videos will be produced.
%
% The annotations on the videos are those described in the documentation of
% 'trackInVideo()' (refer to 'trackInVideo.m').
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
% - R. Hartley and A. Zisserman. Multiple View Geometry in Computer Vision,
%   2nd Edition. Cambridge, UK: Cambridge University Press, 2003.

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created January 19, 2018

%% Input data and parameters

% List of parameters to save with results
parameters_list = {
        'detection_model_filename',...
        'rgb_sigma_filename',...
        'camera_params_filename',...
        'video_filenames',...
        'use_kernel_estimators',...
        'detectionParams',...
        'localizationParams',...
        'options'...
    };

% Probe detection model
detection_model_filename = '/home/llanos/GoogleDrive/PointProbing/DataAndResults/20180114_bluePenWithTape_saturationThreshold/probeDetectionModel_sat0.25.mat';
% RGB noise parameters
rgb_sigma_filename = '/home/llanos/GoogleDrive/PointProbing/DataAndResults/20180112_bluePenWithTape/noiseEstimation/rgbstddev_nonInteractive_video.mat';
% Camera calibration
camera_params_filename = '/home/llanos/GoogleDrive/PointProbing/DataAndResults/20180112_bluePenWithTape/cameraCalibration/cameraParams.mat';

% Wildcard for 'ls()' to find the videos to process.
% Leave empty (`[]`) to read live video
input_video_wildcard = [];

% Output directory for annotated videos
% Leave empty (`[]`) for no output video
output_video_directory = [];

% Output directory for numerical results
% Leave empty (`[]`) for no output data files
output_data_directory = [];

% Video processing options
options.silent = false;

% Parameters which do not usually need to be changed
run('SetFixedParameters.m')

%% Find the videos

if isempty(input_video_wildcard)
    video_filenames = {[]};
    n_videos = 1;
else
    % Find all filenames
    video_filenames = strtrim(strsplit(ls(input_video_wildcard), {'\f','\n','\r','\t','\v'}));
    n_videos = length(video_filenames) - 1; % There is always a terminating newline
end

%% Load calibration data

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

load(camera_params_filename, 'cameraParams');
if ~exist('cameraParams', 'var')
    error('No camera calibration (`cameraParams`) found in ''%s''.', camera_params_filename)
end

%% Process each video (or live video)

save_variables_list = [ parameters_list, {...
        'model_filename',...
        'localizations',...
        'detections'
    } ];

for i = 1:n_videos
    % Output filename for the video
    if (~isempty(output_data_directory) || ~isempty(output_data_directory))...
            && isempty(video_filenames{i})
        cdate = replace(datestr(now, 31), {'-',' ',':'},'_');
    elseif ~isempty(output_data_directory) || ~isempty(output_data_directory)
        [filepath, name, ext] = fileparts(video_filenames{i});
    end
    if isempty(output_video_directory)
        save_video_filename = [];
    else
        if isempty(video_filenames{i})
            save_video_filename = ['live_' cdate '.avi'];
        else
            save_video_filename = [name '_output' ext];
        end
        save_video_filename = fullfile(output_video_directory, save_video_filename);
    end
    
    % Video processing
    [ localizations, detections ] = trackInVideo(...
        video_filenames{i}, save_video_filename,...
        probe, probe_color_distributions, rgb_sigma_polyfit,...
        cameraParams, detectionParams, localizationParams, options...
    );
    
    % Save data to a file
    if ~isempty(output_data_directory)
        if isempty(video_filenames{i})
            save_data_filename = ['live_' cdate '.mat'];
        else
            save_data_filename = [name '_output.mat'];
        end
        save(save_variables_list, fullfile(output_data_directory, save_data_filename));
    end
end
