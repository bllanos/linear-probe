%% Probe Object 3D Localization
% Use the measured dimensions of the probe object and the locations of
% points along its bands in the image to determine the 3D location of the
% probe object.
%
% ## Usage
%   Modify the parameters, and the paths to input data, in the first code
%   section below, then run.
%
% ## Probe Localization Assumptions and Limitations
% - [March 2, 2017] Presently, all interest points on the probe are assumed
%   to be valid, and correctly associated with the measured points on the
%   probe.
%
% ## Input
%
% ### Probe detection result
%
% A '.mat' file containing the variables
% 'probe_detection_matches_filtered', and 'model_filename', output by
% '.\DetectProbe.m'. Refer to the documentation of '.\DetectProbe.m' for
% details.
%
% ### Camera projection parameters
%
% The 3 x 4 camera projection matrix, 'P', provided in a '.mat' file.
%
% ### Image containing the probe
%
% The 'I_filename' variable from the '.mat' file produced by
% '.\DetectProbe.m' is assumed to contain the path to the image in which
% the probe was detected. This image is used to assist with probe
% localization error visualization, but is not necessary for localizing the
% probe.
%
% The image must be an RGB image, in any format that can be loaded by
% `imread`.
%
% ## Output
%
% ### Probe localization results
% A '.mat' file containing the following variables:
%
% - 'model_filename': The path to the '.mat' file containing user-provided
%   measurements of the probe in the structure 'probe'. The
%   'model_filename' parameter is provided to this script via the output of
%   '.\DetectProbe.m', and so 'model_filename' is actually an indirect
%   parameter of this script. It is retrieved from the output of
%   '.\DetectProbe.m' and copied to the output of this script for
%   completeness.
%
% - 'probe_band_locations': A structure vector describing the positions of
%   the detected interest points on the probe. Each element contains the
%   locations of the two endpoints of one detected edge between probe
%   bands. Only probe bands matched to measured positions on the physical
%   probe are included; This structure is based on
%   'probe_detection_matches_filtered' output by '.\DetectProbe.m'. The
%   fields of the structure vector are as follows:
%   - index: The index of the detected edge between two bands of colour
%       on the probe. The indices of edges correspond to the order of the
%       edges along the major axis of the probe in the image, from left to
%       right.
%   - pointAbovePCAMajorAxis: A two-element row vector containing the detected
%       position of the point where the edge meets the lower border of the
%       probe in the image. (This is the point "above" the axis of the
%       probe in the sense that it has a pixel y-coordinate that is larger
%       than that of the axis at the same pixel x-coordinate.)
%   - pointBelowPCAMajorAxis: Similar to 'pointAbovePCAMajorAxis', but contains
%       the detected position of the point where the edge meets the upper
%       border of the probe in the image.
%   - matchedLengthIndex: The index of the edge between coloured bands of
%       the physical probe that is matched with the detected edge.
%       Specifically, this is the index into the user-provided measurements
%       of the probe, `probe.lengths` and `probe.widths`, where `probe` is
%       a structure saved in the file referred to by 'model_filename'.
%
% - 'probe_axis_locations': A structure vector describing the estimated
%   positions of the points along the central axis of the probe object.
%   'probe_band_locations' is useful for assessing reprojection error,
%   whereas 'probe_axis_locations' provides data for further applications
%   (e.g. producing a point cloud by localizing the probe tip in a sequence
%   of probe positions). Each element contains the locations of a point on
%   the line through the centre of the probe, corresponding to a length
%   measurement in `probe.lengths` (`probe` is a structure saved in the
%   file referred to by 'model_filename'). The points are in order, so the
%   first point is the tip of the probe, while the last point is the back
%   end of the probe. The fields of the structure vector are as follows:
%   - imagePoint: A two-element row vector containing the estimated
%       position of the point in the image.
%   - objectPoint: A three-element row vector containing the estimated
%       position of the point in 3D space. (The point is in the coordinate
%       frame determined by the extrinsic parameters of the input camera
%       matrix, 'P'.)
%
% - 'probe_axis': A 3D unit row vector containing the estimated direction
%   of the probe object in space. 'probe_axis' points from the tip of the
%   probe (corresponding to the first measurement in `probe.lengths`) to
%   the end of the probe.
%
% Additionally, the output file contains the values of all parameters in
% the first section of the script below, for reference. (Specifically,
% those listed in `parameters_list`, which should be updated if the set of
% parameters is changed.)
%
% ## References
% - R. Hartley and A. Zisserman. Multiple View Geometry in Computer Vision,
%   2nd Edition. Cambridge, UK: Cambridge University Press, 2003.
%
% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created March 2, 2017

%% Input data and parameters

% List of parameters to save with results
parameters_list = {
        'detection_filename',...
        'camera_filename'...
    };
% The following variables are loaded from the file pointed to by
% 'detection_filename':
parameters_list = [parameters_list {
        'model_filename',...
        'I_filename'...
    }];

% Probe detection result
detection_filename = 'C:\Users\llanos\Google Drive\PointProbing\Data and results\20160811_bambooSkewerProbe\20170217_bambooSkewer_orangeBlue_probeDetectionResults_bottomCamera_rect.mat';

% Camera projection matrix
camera_filename = 'C:\Users\llanos\Google Drive\PointProbing\Data and results\20160811_bambooSkewerProbe\camera_calibration_from_20160609\20170305_bottomCameraMatrix_identityExtrinsics.mat';

%% Load input data
detection_variables_required = {...
        'probe_detection_matches_filtered',...
        'model_filename',...
        'I_filename'...
    };
load(detection_model_filename, detection_variables_required{:});
if ~all(ismember(detection_variables_required, who))
    error('One or more of the probe detection variables is not loaded.')
end

load(camera_filename, 'P');
if ~exist('P', 'var')
    error('No camera matrix found in ''%s''.',camera_filename)
end

I = imread(I_filename);