%% Evaluation of Probe Object 3D Localization
% Assess the estimated 3D pose of the probe object in terms of its
% reprojection error in an image.
%
% ## Usage
%   Modify the parameters, and the paths to input data, in the first code
%   section below, then run.
%
% ## Input
%
% ### Probe localization result
%
% A '.mat' file containing the variable 'probe_axis_locations', and,
% optionally, the following parameters (which can be overridden by defining
% them in this script):
% - model_filename
% - camera_filename
% - I_filename
%
% These items have the same form as output by '.\LocalizeProbe.m'. Refer to
% the documentation of '.\LocalizeProbe.m' for details.
%
% The path to this '.mat' file is provided by the 'localization_filename'
% parameter in this script.
%
% ### Probe measurements
%
% The 'model_filename' parameter is the path to a '.mat' file containing a
% structure called 'probe'. 'probe' is a set of user-provided measurements
% of the physical probe. Refer to the documentation of
% '.\CreateProbeDetectionModel.m' for details.
%
% ### Camera projection parameters
%
% The 3 x 4 camera projection matrix, 'P', provided in a '.mat' file. The
% path to the '.mat' file is provided in the 'camera_filename' parameter.
%
% Note: The camera should be calibrated using the same units of measurement
% as used for the measurements of the probe in 'model_filename'.
%
% ### Image containing the probe
%
% The 'I_filename' parameter is assumed to contain the path to the image
% showing the probe in the location corresponding to localization result.
% The image should correspond to the camera 'P', but this camera may not
% have provided the image in which the probe was detected.
%
% The image must be an RGB image, in any format that can be loaded by
% `imread`.
%
% ### Annotations for probe image (Optional)
%
% If this data, the ground truth, is provided, localization error can be
% quantified.
%
% Data is provided in a '.mat' file whose filepath is stored in the
% parameter 'visible_points_filename'.
%
% The '.mat' file contains the variable 'visible_points', which is an n x 3
% array, where the columns represent the following:
% - visible_points(:, 1): Indices in 'probe_axis_locations' corresponding
%     to the user-marked interest point.
% - visible_points(:, 2): Image x-coordinates of user-marked interest
%     points
% - visible_points(:, 3): Image y-coordinates of user-marked interest
%     points
%
% The user-marked points are the ground truth locations on the contour of
% the probe corresponding to the points in 'probe_axis_locations'. (To be
% specific, the user-marked points are the points on the probe contour
% on the ground truth cross-sections corresponding to the points in
% 'probe_axis_locations'.)
%
% A single interest point should be marked for an end of the probe that
% tapers to a point. Two interest points should be marked for
% cross-sections of the probe that have non-negligible widths. The system
% will automatically determine which points are the 'bottom' and 'top'
% points when assessing reprojection error.
%
% ## Output
%
% No output will be produced in the absence of ground truth input data.
%
% ### Probe localization evaluation results
% A '.mat' file containing the following variables:
%
% - 'visible_points_error': A structure vector comparing the reprojected
%   and user-marked points. The fields of the structure vector are as
%   follows:
%   - index: The index of the cross-section of the probe in
%       'probe_axis_locations'. The elements of the structure vector are
%       sorted by this field.
%   - gtAboveAxis: A two-element row vector containing the user-marked
%       position of the interest point on the lower border of the
%       probe in the image. (This is the point "above" the axis of the
%       probe in the sense that it has a pixel y-coordinate that is larger
%       than that of the axis at the same pixel x-coordinate.)
%   - gtBelowAxis: Similar to 'gtAboveAxis', but contains the user-marked
%       position of the interest point on the upper border of the
%       probe in the image.
%   - reprojectedAboveAxis: A two-element row vector containing the
%       reprojected position of the interest point on the lower border of
%       the probe in the image.
%   - reprojectedBelowAxis: Similar to 'reprojectedAboveAxis', but contains
%       the reprojected position of the interest point on the upper border
%       of the probe in the image.
%   - error: The sum of the reprojection errors (Euclidean distances
%       between reprojected and ground-truth points) for the two points on
%       the cross-section of the probe.
%   Ground truth points for the tapered tip(s) of the probe are simply
%   duplicated so that both the 'gtAboveAxis' and 'gtBelowAxis' fields are
%   populated for all locations along the probe. However, the 'error' field
%   for a tapered tip is computed from a single distance measurement.
%
% - 'mean_error': The average of the 'error' field values in
%   'visible_points_error' (avoiding double-counting the probe tips).
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
% File created March 27, 2017

%% Input data and parameters

% List of parameters to save with results
parameters_list = {
        'localization_filename',...
        'visible_points_filename'...
    };
% The following variables are loaded from the file pointed to by
% 'localization_filename', unless they are defined below:
parameters_list = [parameters_list {
        'model_filename',...
        'camera_filename',...
        'I_filename'...
    }];

% Probe localization results
localization_filename = 'C:\Users\llanos\Google Drive\PointProbing\Data and results\20160811_bambooSkewerProbe\20170409_bambooSkewer_orangeBlue_probeLocalizationResults_bottomCamera_rect.mat';

% Measurements of the probe object
model_filename = [];

% Camera projection matrix
camera_filename =  [];

% Image containing the probe
I_filename = [];

% Ground truth points
% If none is provided, set it to an empty array (`[]`).
visible_points_filename = [];

% Flags to control verbose or graphical output
display_reprojection = true;

%% Load input data
localization_variables_required = {...
        'probe_axis_locations'...
    };
if ~exist('model_filename', 'var') || isempty(model_filename)
    localization_variables_required = [localization_variables_required 'model_filename'];
end
if ~exist('camera_filename', 'var') || isempty(I_filename)
    localization_variables_required = [localization_variables_required 'camera_filename'];
end
if ~exist('I_filename', 'var') || isempty(I_filename)
    localization_variables_required = [localization_variables_required 'I_filename'];
end
load(localization_filename, localization_variables_required{:});
if ~all(ismember(localization_variables_required, who))
    error('One or more of the probe localization variables is not loaded.')
end

load(model_filename, 'probe');
if ~exist('probe', 'var')
    error('No variable called ''probe'' loaded from %s (which would contain probe measurements).', model_filename)
end

load(camera_filename, 'P');
if ~exist('P', 'var')
    error('No camera matrix found in ''%s''.',camera_filename)
end

I = imread(I_filename);
image_size = size(I);
image_size = image_size(1:2);

single_view_gt_available = exist('visible_points_filename', 'var') && ~isempty(visible_points_filename);
if single_view_gt_available
    load(visible_points_filename, 'visible_points');
    if ~exist('visible_points', 'var')
        error('No variable called ''visible_points'' loaded from ''%s''.',visible_points_filename)
    end
end

%% Reproject the probe into the image

X_tip = probe_axis_locations(1).objectPoint;
X_end = probe_axis_locations(end).objectPoint;
[d, ~] = probeOrientationFromEndpoints( P, X_tip, X_end );

if display_reprojection
    if single_view_gt_available
        plotProbeReprojection(...
                    I, visible_points(:, 2:3), probe.lengths, probe.widths, P, d, X_tip,...
                    'Reprojection of probe localization result'...
                );
    else
        plotProbeReprojection(...
            I, [], probe.lengths, probe.widths, P, d, X_tip,...
            'Reprojection of probe localization result'...
        );
    end
end