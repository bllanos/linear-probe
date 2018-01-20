function [ localizations, detections ] = trackInVideo(...
    in_filename, out_filename,...
    probe, probe_color_distributions, rgb_sigma_polyfit,...
    cameraParams, detectionParams, localizationParams, options...
    )
% TRACKINVIDEO  Track the probe in video data
%
% ## Syntax
% localizations = trackInVideo(...
%     in_filename, out_filename,...
%     probe, probe_color_distributions, rgb_sigma_polyfit,...
%     cameraParams, detectionParams, localizationParams, options...
% )
% [ localizations, detections ] = trackInVideo(...
%     in_filename, out_filename,...
%     probe, probe_color_distributions, rgb_sigma_polyfit,...
%     cameraParams, detectionParams, localizationParams, options...
% )
%
% ## Description
% localizations = trackInVideo(...
%     in_filename, out_filename,...
%     probe, probe_color_distributions, rgb_sigma_polyfit,...
%     cameraParams, detectionParams, localizationParams, options...
% )
%   Return probe localization results in video frames.
% [ localizations, detections ] = trackInVideo(...
%     in_filename, out_filename,...
%     probe, probe_color_distributions, rgb_sigma_polyfit,...
%     cameraParams, detectionParams, localizationParams, options...
% )
%   Additionally return probe detection results in video frames.
%
% ## Input Arguments
%
% in_filename -- Input video
%   A character vector containing the name and path of the input video
%   file. If empty (`[]`), the function will attempt to read live video
%   from a webcam, using the MATLAB Support Package for USB Webcams.
%
% out_filename -- Output video
%   A character vector containing the name and path of the output video
%   file, which is an annotated version of the input video. All frames
%   which are in the input video or which were captured from live video
%   will be output. Frames in which the probe was successfully localized
%   will be annotated by 'plotProbeReprojection()'. For detection results
%   which are not sufficient for localizing the probe, only the detected
%   points will be annotated.
%
%   If empty (`[]`), no output video will be produced.
%
% probe -- Probe measurements
%   Refer to the documentation of 'CreateProbeDetectionModel.m' for
%   details.
%
% probe_color_distributions -- Probe colour estimators
%   Discretized density estimators of image hue values corresponding to the
%   different coloured bands on the probe, in the same order (starting from
%   the active tip of the probe). The i-th column of this 2D array stores
%   the estimator for the i-th colour class of probe segments.
%
% rgb_sigma_polyfit -- Camera RGB noise model
%   A parameter passed indirectly to 'extractBinaryRegions()'. As described
%   in the documentation of 'extractBinaryRegions()'. An empty array can be
%   passed instead, if `detectionParams.uniform_background_initial` and
%   `detectionParams.uniform_background_final` are both `false`.
%
% cameraParams -- Camera calibration parameters
%   A structure of class 'cameraParameters' describing the camera, such as
%   output by the MATLAB Camera Calibrator app. Only the 'IntrinsicMatrix',
%   'RadialDistortion', and 'TangentialDistortion' fields are used. The
%   'IntrinsicMatrix' field is used to construct the camera matrix with the
%   camera serving as the origin of the coordinate frame. The
%   'RadialDistortion', and 'TangentialDistortion' fields are implicitly
%   used to undistort image coordinates by MATLAB's 'undistortPoints()'
%   and 'undistortImage()' function.
%
% detectionParams -- Probe detection fixed parameters
%   The `params` input argument of 'detectProbe()'.
%
% localizationParams -- Probe localization fixed parameters
%   The `params` input argument of 'localizeProbe()'.
%
% options -- Processing options
%   A structure with the following field:
%   - silent: If `true`, a video player will not be opened to show the
%     video frames annotated with probe detection and localization results.
%
% ## Output Arguments
%
% localizations -- Probe localization results
%   A structure vector, where each element describes the localization of
%   the probe in a video frame. There are no entries for video frames in
%   which the probe was not localized. Each element has the following
%   fields:
%   - frame: The index of the frame (starting from `1`) corresponding to
%     the localization result.
%   - axis_locations: The `axis_locations` output argument of
%     'localizeProbe()'
%   - probe_axis: The `probe_axis` output argument of
%     'localizeProbe()'
%   - band_locations: The `band_locations` output argument of
%     'localizeProbe()'
%
% detections -- Probe detection results
%   A structure vector, where each element describes the detection of the
%   probe in a video frame. There are no entries for video frames in which
%   the probe was not detected. Note that there may be some frames for
%   which there are elements in `detections`, but no elements in
%   `localizations`. Each element has the following fields:
%   - frame: The index of the frame (starting from `1`) corresponding to
%     the detection result.
%   - matches_filtered: The `matches_filtered` output argument of
%     'detectProbe()'
%   - matches: The `matches` output argument of
%     'detectProbe()'
%
% ## Notes
% - Image points are undistorted, rather than entire images, prior to probe
%   localization. The amount of distortion is assumed to be sufficiently
%   small not to interfere with the function of 'detectProbe()'. However,
%   in order to produce correct annotated output images, entire images are
%   undistorted prior to annotation.
% - All image coordinates in `detections` and `localizations` are corrected
%   for lens distortion.
%
% ## References
% - MATLAB Example: Face Detection and Tracking Using Live Video
%   Acquisition
%
% See also detectProbe, localizeProbe, VideoReader, undistortPoints, plotProbeReprojection

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created January 19, 2018

nargoutchk(1,2)
narginchk(9,9)

use_live_video = isempty(in_filename);
if use_live_video
    inputVideo = webcam();
else
    inputVideo = vision.VideoFileReader(in_filename);
end

if use_live_video
    I = snapshot(inputVideo);
else
    I = step(inputVideo);
end
image_size = size(I);

if ~options.silent
    player = vision.VideoPlayer('Position', [100 100 [image_size(2), image_size(1)]+30]);
elseif use_live_video
    error('There is no way to gracefully stop live video capture without opening a video player.')
end

runLoop = true;
while runLoop
    if ~options.silent
        step(player, I);
    end
    
    if ~use_live_video
        runLoop = ~isDone(inputVideo);
    end
    if ~options.silent
        runLoop = runLoop && isOpen(player);
    end
    if runLoop
        if use_live_video
            I = snapshot(inputVideo);
        else
            I = step(inputVideo);
        end
    end
end

% Cleanup
if use_live_video
    clear inputVideo;
else
    release(inputVideo);
end
if ~options.silent
    release(player);
end

localizations = [];
detections = [];

end