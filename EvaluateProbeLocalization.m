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
% optionally, 'model_filename'. These variables have the same form as
% output by '.\LocalizeProbe.m'. Refer to the documentation of
% '.\LocalizeProbe.m' for details.
%
% ### Probe measurements
%
% The 'model_filename' parameter is the path to a '.mat' file containing a
% structure called 'probe'. 'probe' is a set of user-provided measurements
% of the physical probe. Refer to the documentation of
% '.\CreateProbeDetectionModel.m' for details.
%
% This parameter is to be provided separately if not contained in the probe
% localization results described above.
%
% ### Camera projection parameters
%
% The 3 x 4 camera projection matrix, 'P', provided in a '.mat' file.
%
% Note: The camera should be calibrated using the same units of measurement
% as used for the measurements of the probe in 'model_filename'.
%
% ### Image containing the probe
%
% The 'I_filename' variable is assumed to contain the path to the image
% showing the probe in the location corresponding to localization result.
% The image should correspond to the camera 'P', but this camera may not
% have provided the image in which the probe was detected.
%
% The image must be an RGB image, in any format that can be loaded by
% `imread`.
%
% ### Annotations for probe image
%
% A '.mat' file containing the variable 'visible_points', which is an n x 3
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
% ### Probe localization evaluation results
% A '.mat' file containing the following variables:
%
% - 'model_filename': The path to the '.mat' file containing user-provided
%   measurements of the probe in the structure 'probe'. The
%   'model_filename' parameter may be provided to this script via the
%   output of '.\LocalizeProbe.m' (see above), in which case
%   'model_filename' is actually an indirect parameter of this script. It
%   is copied to the output of this script for completeness.
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
%   'visible_points_error'.
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