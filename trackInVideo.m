function [ localizations, detections ] = trackInVideo(...
    in_filename, out_filenames,...
    probe, probe_color_distributions, rgb_sigma_polyfit,...
    cameraParams, detectionParams, localizationParams, options...
    )
% TRACKINVIDEO  Track the probe in video data
%
% ## Syntax
% trackInVideo(...
%     in_filename, out_filenames,...
%     probe, probe_color_distributions, rgb_sigma_polyfit,...
%     cameraParams, detectionParams, localizationParams, options...
% )
% localizations = trackInVideo(...
%     in_filename, out_filenames,...
%     probe, probe_color_distributions, rgb_sigma_polyfit,...
%     cameraParams, detectionParams, localizationParams, options...
% )
% [ localizations, detections ] = trackInVideo(...
%     in_filename, out_filenames,...
%     probe, probe_color_distributions, rgb_sigma_polyfit,...
%     cameraParams, detectionParams, localizationParams, options...
% )
%
% ## Description
% trackInVideo(...
%     in_filename, out_filenames,...
%     probe, probe_color_distributions, rgb_sigma_polyfit,...
%     cameraParams, detectionParams, localizationParams, options...
% )
%   Track the probe in video frames, and output results to files.
% localizations = trackInVideo(...
%     in_filename, out_filenames,...
%     probe, probe_color_distributions, rgb_sigma_polyfit,...
%     cameraParams, detectionParams, localizationParams, options...
% )
%   Additionally return probe localization results in video frames.
% [ localizations, detections ] = trackInVideo(...
%     in_filename, out_filenames,...
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
%   from either a webcam, using the MATLAB Support Package for USB Webcams,
%   or a Point Grey camera, using the Image Acquisition Toolbox Support
%   Package for Point Grey Hardware. The choice of camera is determined by
%   the `options.video_mode` parameter (see below).
%
% out_filenames -- Output filepaths
%   A structure with the following fields. Each field is a character vector
%   containing a file name and path. If a field is empty (`[]`), the
%   corresponding file will not be produced:
%   - raw_video: A copy of the input video (not useful except for live
%     video capture).
%   - out_video: An annotated version of the input video. All frames
%     which are in the input video, or which were captured from live video
%     will be output. Frames in which the probe was successfully localized
%     will be annotated by 'plotProbeReprojection()'. For detection results
%     which are not sufficient for localizing the probe, only the detected
%     points will be annotated. All frames output to this video are
%     corrected for lens distortion.
%   - out_csv: A CSV file with the following columns. One row is generated
%     for each frame in which the probe was successfully localized.
%     - Video index (the value of `options.video_index`)
%     - Frame index (the first frame has index `1`)
%     - x-coordinate of the tip of the probe in 3D space
%     - y-coordinate of the tip of the probe in 3D space
%     - z-coordinate of the tip of the probe in 3D space
%     - x-coordinate of the unit direction vector from the tip to the end
%       of the probe in 3D space
%     - y-coordinate of the unit direction vector from the tip to the end
%       of the probe in 3D space
%     - z-coordinate of the unit direction vector from the tip to the end
%       of the probe in 3D space
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
%   Extrinsic parameters in the structure are ignored because they most
%   likely originate from irrelevant camera calibration images.
%
% detectionParams -- Probe detection fixed parameters
%   The `params` input argument of 'detectProbe()'.
%
% localizationParams -- Probe localization fixed parameters
%   The `params` input argument of 'localizeProbe()'.
%
% options -- Processing options
%   A structure with the following fields:
%   - video_mode: When video input is loaded from files, this field is
%     ignored. In the case of live video input, if this field is 'webcam',
%     the function will attempt to connect to a USB webcam. Otherwise, it
%     will use the value of this field as the `format` parameter of
%     'videoinput()', and will attempt to connect to a Point Grey camera.
%   - silent: If `true`, a video player will not be opened to show the
%     video frames annotated with probe detection and localization results.
%   - frame_rate: The output videos' framerates.
%   - record_only: If `true`, just copy the video to the output file
%     without processing. In this case, `localizations` and `detections`
%     will be empty, and the only output file that can be generated is a
%     copy of the input video. The video player will also display raw
%     images, rather than undistorted images. For live video, more frames
%     can be captured per unit time than if processing was enabled.
%   - show_errors: Output exceptions, which are thrown by the probe
%     detection and localization functions, as warnings. These exceptions
%     are bugs - Exceptions thrown because of failed detection are not
%     output, for instance.
%   - video_index: The value to be output in the first column for each row
%     of any output CSV file. This field is not required if
%     `out_filenames.out_csv` is empty.
%   - append_csv: A Boolean value indicating whether to overwrite the
%     output CSV file (`false`) or append to it (`true`). This field is not
%     required if `out_filenames.out_csv` is empty.
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
%   - hyp: The `hyp` output argument of
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
%     'detectProbe()', with image point coordinates corrected for lens
%     distortion.
%   - matches: The `matches` output argument of
%     'detectProbe()'. Image point coordinates are not corrected for lens
%     distortion.
%
% ## Notes
% - Image points are undistorted, rather than entire images, prior to probe
%   localization. The amount of distortion is assumed to be sufficiently
%   small not to interfere with the function of 'detectProbe()'. However,
%   in order to produce correct annotated output images, entire images are
%   undistorted prior to annotation.
% - All image coordinates in `localizations` are corrected for lens
%   distortion.
%   
% ## References
% - MATLAB Example: Face Detection and Tracking Using Live Video
%   Acquisition
%
% See also detectProbe, localizeProbe, VideoReader, undistortPoints, plotProbeReprojection, webcam, videoinput

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created January 19, 2018

nargoutchk(0,2)
narginchk(9,9)

% Connect to video source
use_live_video = isempty(in_filename);
if use_live_video
    use_webcam = strcmp(options.video_mode, 'webcam');
    if use_webcam
        inputVideo = webcam();
    else
        inputVideo = videoinput('pointgrey', 1, options.video_mode);
    end
else
    inputVideo = vision.VideoFileReader(in_filename);
end
if use_live_video
    if use_webcam
        I = snapshot(inputVideo);
    else
        I = getsnapshot(inputVideo);
    end
else
    I = step(inputVideo);
end
frame_index = 1;
past_frame_index = 0;
console_output_interval = 3;

% Initialize video player
image_size = size(I);
image_size = image_size(1:2);
[max_image_size, max_image_size_ind] = max(image_size);
player_max = 600;
if max_image_size > player_max
    player_scale = player_max / image_size(max_image_size_ind);
else
    player_scale = 1;
end
player_size = ceil(image_size * player_scale);
if ~options.silent
    player = vision.VideoPlayer(...
        'Position',...
        [100 100 (fliplr(player_size) +30)]...
    );
elseif use_live_video
    error('There is no way to gracefully stop live video capture without opening a video player.')
end

% Initialize output streams
raw_video_output_enabled = ~isempty(out_filenames.raw_video);
if raw_video_output_enabled
    outputRawVideo = vision.VideoFileWriter(out_filenames.raw_video,'FrameRate',options.frame_rate);
end
if ~options.record_only
    annotated_video_output_enabled = ~isempty(out_filenames.out_video);
    annotate_video = annotated_video_output_enabled || ~options.silent; 
    if annotated_video_output_enabled
        outputAnnotatedVideo = vision.VideoFileWriter(out_filenames.out_video,'FrameRate',options.frame_rate);
    end
    csv_output_enabled = ~isempty(out_filenames.out_csv);
    if csv_output_enabled
        if options.append_csv
            outputCSV = fopen(out_filenames.out_csv, 'a');
        else
            outputCSV = fopen(out_filenames.out_csv, 'w');
        end
    end
end

% Process video
runLoop = true;
if ~options.record_only
    return_localizations = (nargout > 0);
    if return_localizations
        localizations = struct(...
            'frame', {}, 'axis_locations', {}, 'probe_axis', {},...
            'band_locations', {}, 'hyp', {}...
        );
    end
    return_detections = (nargout > 1);
    if return_detections
        detections = struct(...
            'frame', {}, 'matches_filtered', {}, 'matches', {}...
        );
    end

    P = [cameraParams.IntrinsicMatrix.' zeros(3, 1)];
    if annotate_video
        marker_size = ceil(max(5, 0.01 * max(image_size)));
    end
end
tic;
while runLoop
    if ~options.record_only
        % Undistort image
        if annotate_video
            I_out = undistortImage(I, cameraParams);
        else
            I_out = I;
        end
            
        try
            % Detect in current frame
            [ matches_filtered, matches ] = detectProbe(...
                   I, probe, probe_color_distributions, rgb_sigma_polyfit, detectionParams...
            );

            % Undistort detected point coordinates
            n_hypotheses = length(matches_filtered);
            any_detected = false;
            for hyp = 1:n_hypotheses
                above = vertcat(matches_filtered{hyp}(:).pointAbovePCAMajorAxis);
                below = vertcat(matches_filtered{hyp}(:).pointBelowPCAMajorAxis);
                n_points = size(above, 1);
                if n_points > 0
                    any_detected = true;

                    undistortedPoints = undistortPoints([ above; below ], cameraParams);
                    above = num2cell(undistortedPoints(1:n_points, :), 2);
                    below = num2cell(undistortedPoints((n_points + 1):end, :), 2);
                    [matches_filtered{hyp}(:).pointAbovePCAMajorAxis] = above{:};
                    [matches_filtered{hyp}(:).pointBelowPCAMajorAxis] = below{:};
                end
            end

            if return_detections && any_detected
                detections(end + 1).matches_filtered = matches_filtered; %#ok<AGROW>
                detections(end).matches = matches;
                detections(end).frame = frame_index;
            end

            if n_points > 2
                % Localize in current frame
                [ axis_locations, probe_axis, band_locations, hyp ] = localizeProbe(...
                    probe, matches_filtered, P, localizationParams...
                );
                if return_localizations
                    localizations(end + 1).axis_locations = axis_locations; %#ok<AGROW>
                    localizations(end).probe_axis = probe_axis;
                    localizations(end).band_locations = band_locations;
                    localizations(end).hyp = hyp;
                    localizations(end).frame = frame_index;
                end

                tip = axis_locations(1).objectPoint;

                if csv_output_enabled
                    fprintf(...
                        outputCSV, '%d, %d, %g, %g, %g, %g, %g, %g\n',...
                        options.video_index, frame_index, tip(1), tip(2), tip(3),...
                        probe_axis(1), probe_axis(2), probe_axis(3)...
                    );
                end

                if annotate_video
                    lengths = vertcat(matches_filtered{hyp}(:).matchedLength);
                    widths = vertcat(matches_filtered{hyp}(:).matchedWidth);
                    I_out = plotProbeReprojection(...
                        I_out, undistortedPoints, lengths, widths, P, probe_axis, tip...
                    );
                end

            elseif n_points > 0 && annotate_video
                points_clipped = undistortedPoints(...
                    undistortedPoints(:, 1) >= 1 & undistortedPoints(:, 2) >= 1 &...
                    undistortedPoints(:, 1) <= image_size(2) & undistortedPoints(:, 2) < image_size(1),...
                    :...
                );
                I_out = insertMarker(...
                    I_out, points_clipped, 'o', 'Color', 'green', 'size', marker_size...
                );
            end
            
        catch ME
            switch ME.identifier
                case 'PROBE:InsufficientRegionsDetected'
                otherwise
                    % rethrow(ME)
                    if options.show_errors
                        warning('Error in video %d, frame %d: "%s"', options.video_index, frame_index, getReport(ME));
                    end
            end
        end
        
        if annotated_video_output_enabled
            step(outputAnnotatedVideo, I_out);
        end
    else
        I_out = I;
    end
    
    if ~options.silent
        step(player, imresize(I_out, player_scale));
    end
    
    if raw_video_output_enabled
        step(outputRawVideo, I);
    end
    
    % Check whether to end processing
    if ~use_live_video
        runLoop = ~isDone(inputVideo);
    end
    if ~options.silent
        runLoop = runLoop && isOpen(player);
    end
    if runLoop
        if use_live_video
            if use_webcam
                I = snapshot(inputVideo);
            else
                I = getsnapshot(inputVideo);
            end
        else
            I = step(inputVideo);
        end
        
        elapsed = toc;
        if elapsed >= console_output_interval
            fprintf(...
                '%d frames processed, %g FPS for the last %g seconds.\n',...
                frame_index, (frame_index - past_frame_index) / elapsed, elapsed...
            );
            past_frame_index = frame_index;
            tic;
        end
        frame_index = frame_index + 1;
    end
end

% Cleanup
if use_live_video
    if use_webcam
        clear inputVideo;
    else
        delete(inputVideo);
    end
else
    release(inputVideo);
end
if ~options.silent
    release(player);
end
if raw_video_output_enabled
    release(outputRawVideo);
end
if ~options.record_only
    if annotated_video_output_enabled
        release(outputAnnotatedVideo);
    end
    if csv_output_enabled
        fclose(outputCSV);
    end
end

if options.record_only
    localizations = [];
    detections = [];
end

end