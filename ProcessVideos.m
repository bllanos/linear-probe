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
% One output file of each of the following types will be generated for each
% input video source. Output files will be saved in the directories
% referred to by the 'output_*_directory' variables below. If an output
% directory variable is empty, no output data files of the corresponding
% type will be produced.
%
% ### Input videos
%
% Copies of the input videos, or videos captured from cameras, will be
% saved in the directory referred to by 'output_raw_video_directory'
% below. These videos are not corrected for lens distortion.
%
% ### Annotated videos
%
% The input videos will be annotated with the reprojection of the probe,
% and saved in the directory referred to by
% 'output_annotated_video_directory' below. These videos are corrected for
% lens distortion.
%
% The annotations on the videos are those described in the documentation of
% 'trackInVideo()' (refer to 'trackInVideo.m').
%
% ### Point clouds
%
% A CSV file will be generated containing frame numbers, the corresponding
% 3D locations of the probe tip, and the corresponding 3D unit direction
% vectors from the probe tip to its other end. The CSV file will be saved
% in the directory referred to by 'output_point_cloud_directory' below.
%
% The output point cloud is also described in the documentation of the
% `out_filenames.out_csv` input argument of 'trackInVideo()' (refer to
% 'trackInVideo.m').
%
% ### Probe detection and localization results
%
% A '.mat' file containing the following variables will be saved in the
% directory referred to by 'output_data_directory' below:
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
% Additionally, the file contains the values of all parameters in the first
% section of the script below, for reference. (Specifically, those listed
% in `parameters_list`, which should be updated if the set of parameters is
% changed.)
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
detection_model_filename = '/home/llanos/GoogleDrive/PointProbing/DataAndResults/20180121_bluePenGreenTape/model/bluePenGreenTape_detectionModel.mat';
% RGB noise parameters
rgb_sigma_filename = '/home/llanos/GoogleDrive/PointProbing/DataAndResults/20180121_bluePenGreenTape/model/rgbstddev_nonInteractive_video.mat';
% Camera calibration
camera_params_filename = '/home/llanos/GoogleDrive/PointProbing/DataAndResults/20180112_bluePenWithTape/cameraCalibration/cameraParams.mat';

% Wildcard for 'ls()' to find the videos to process.
% Leave empty (`[]`) to read live video
input_video_wildcard = [];

% Output directory for raw videos
% Leave empty (`[]`) for no output raw video
output_raw_video_directory = '/home/llanos/Downloads/output';

% Output directory for annotated videos
% Leave empty (`[]`) for no output annotated video
output_annotated_video_directory = '/home/llanos/Downloads/output';

% Output directory for CSV format point cloud
% Leave empty (`[]`) for no output point cloud
output_point_cloud_directory = '/home/llanos/Downloads/output';

% Output directory for comprehensive numerical results
% Leave empty (`[]`) for no output data file
output_data_directory = [];

% Video processing options
% Refer to the documentation of the `options` parameter of 'trackInVideo()'
% in 'trackInVideo.m'.
options.silent = false;
options.frame_rate = 20;
options.record_only = false;
options.show_errors = true;

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
any_output_files = ~all([
        isempty(output_raw_video_directory);
        isempty(output_annotated_video_directory);
        isempty(output_point_cloud_directory);
        isempty(output_data_directory)
    ]);

for i = 1:n_videos
    % Generate output filenames
    if any_output_files && isempty(video_filenames{i})
        cdate = replace(datestr(now, 31), {'-',' ',':'},'_');
    elseif any_output_files
        [filepath, name, ext] = fileparts(video_filenames{i});
    end
    if isempty(output_raw_video_directory)
        raw_video_filename = [];
    else
        if isempty(video_filenames{i})
            raw_video_filename = ['live_' cdate '.avi'];
        else
            raw_video_filename = [name '_copy' ext];
        end
        raw_video_filename = fullfile(output_raw_video_directory, raw_video_filename);
    end
    if isempty(output_annotated_video_directory)
        annotated_video_filename = [];
    else
        if isempty(video_filenames{i})
            annotated_video_filename = ['live_' cdate '_annotated.avi'];
        else
            annotated_video_filename = [name '_annotated' ext];
        end
        annotated_video_filename = fullfile(output_annotated_video_directory, annotated_video_filename);
    end
    if isempty(output_point_cloud_directory)
        point_cloud_filename = [];
    else
        if isempty(video_filenames{i})
            point_cloud_filename = ['live_' cdate '_pointCloud.csv'];
        else
            point_cloud_filename = [name '_pointCloud.csv'];
        end
        point_cloud_filename = fullfile(output_point_cloud_directory, point_cloud_filename);
    end
    if ~isempty(output_data_directory)
        if isempty(video_filenames{i})
            save_data_filename = ['live_' cdate '_fullResults.mat'];
        else
            save_data_filename = [name '_fullResults.mat'];
        end
        save_data_filename = fullfile(output_data_directory, save_data_filename);
    end
    
    out_filenames = struct(...
        'raw_video', raw_video_filename,...
        'out_video', annotated_video_filename,...
        'out_csv', point_cloud_filename...
    );
    
    % Video processing    
    if ~isempty(output_data_directory)
        [ localizations, detections ] = trackInVideo(...
            video_filenames{i}, out_filenames,...
            probe, probe_color_distributions, rgb_sigma_polyfit,...
            cameraParams, detectionParams, localizationParams, options...
        );
        save(save_data_filename, save_variables_list{:});
    else
        trackInVideo(...
            video_filenames{i}, out_filenames,...
            probe, probe_color_distributions, rgb_sigma_polyfit,...
            cameraParams, detectionParams, localizationParams, options...
        );
    end
end
