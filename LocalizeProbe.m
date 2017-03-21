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
% 'probe_detection_matches_filtered', and, optionally, 'model_filename',
% output by '.\DetectProbe.m'. Refer to the documentation of
% '.\DetectProbe.m' for details.
%
% ### Probe measurements
%
% The 'model_filename' parameter is the path to a '.mat' file containing a
% structure called 'probe'. 'probe' is a set of user-provided measurements
% of the physical probe. Refer to the documentation of
% '.\CreateProbeDetectionModel.m' for details.
%
% This parameter is to be provided separately if not contained in the probe
% detection results described above.
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
% The 'I_filename' variable, from the '.mat' file produced by
% '.\DetectProbe.m', or provided separately, is assumed to contain the path
% to the image in which the probe was detected. This image is used to
% assist with probe localization error visualization, but is not necessary
% for localizing the probe.
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
%   'probe_detection_matches_filtered', output by '.\DetectProbe.m', except
%   that the points are reprojected from their estimated 3D locations. The
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
        'camera_filename',...
        'linear_convergence_threshold',...
        'normalize_homography1D'...
    };
% The following variables are loaded from the file pointed to by
% 'detection_filename', unless they are defined below:
parameters_list = [parameters_list {
        'model_filename',...
        'I_filename'...
    }];

% Measurements of the probe object
model_filename = 'C:\Users\llanos\Google Drive\PointProbing\Data and results\20160811_bambooSkewerProbe\bambooSkewer_orangeBlue.mat';

% Probe detection result
detection_filename = 'C:\Users\llanos\Google Drive\PointProbing\Data and results\20160811_bambooSkewerProbe\20170217_bambooSkewer_orangeBlue_probeDetectionResults_bottomCamera_rect.mat';

% Camera projection matrix
camera_filename = 'C:\Users\llanos\Google Drive\PointProbing\Data and results\20160811_bambooSkewerProbe\camera_calibration_from_20160609\20170305_bottomCameraMatrix_identityExtrinsics.mat';

% Linear Probe Localization
% Error convergence threshold for linear probe estimation
linear_convergence_threshold = 0.01;
% Normalize lengths when estimating a 1D homography between the probe and
% its image
normalize_homography1D = true;

% Debugging tools
verbose_linear_estimation = false; % Requires `I_filename` to be valid
display_linear_estimation = true; % Requires `I_filename` to be valid
verbose_nonlinear_estimation = true;
display_nonlinear_estimation = true; % Requires `I_filename` to be valid

%% Load input data
detection_variables_required = {...
        'probe_detection_matches_filtered'...
    };
if ~exist('model_filename', 'var')
    detection_variables_required = [detection_variables_required 'model_filename'];
end
if ~exist('I_filename', 'var')
    detection_variables_required = [detection_variables_required 'I_filename'];
end
load(detection_filename, detection_variables_required{:});
if ~all(ismember(detection_variables_required, who))
    error('One or more of the probe detection variables is not loaded.')
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

%% Linear estimation of probe location

% Estimated centerline of the probe in the image
above = vertcat(probe_detection_matches_filtered(:).pointAbovePCAMajorAxis);
below = vertcat(probe_detection_matches_filtered(:).pointBelowPCAMajorAxis);
allPoints = [ above; below ];
[ coeff, ~, ~, ~, ~, mu ] = pca(allPoints);

% Express the PCA component vectors as lines in the space of the original data
dx = coeff(1, :).';
dy = coeff(2, :).';
c = ones(2, 1);
a = c ./ ((dx * mu(2) ./ dy) - mu(1));
b = -a .* dx ./ dy;
% Normalize
scale = sqrt(a .^ 2 + b .^ 2);
axes = [a, b, c] ./ repmat(scale, 1, 3);

% The first axis is the estimated centerline of the probe
%
% I have not found a way to estimate the true axis of the probe (taking
% projective distortion into account) that remains numerically stable as
% the width of the probe decreases. However, as the width of the probe
% decreases, the first PCA axis becomes an increasingly good approximation
% of the true probe axis.
image_centerline = axes(1, :);

lengths = vertcat(probe_detection_matches_filtered(:).matchedLength);
widths = vertcat(probe_detection_matches_filtered(:).matchedWidth);
if verbose_linear_estimation
    [X_tip, d, u] = probeTipAndOrientation(...
        above, below, lengths, widths, P, image_centerline,...
        linear_convergence_threshold, normalize_homography1D, I...
    ); %#ok<UNRCH>
else
    [X_tip, d, u] = probeTipAndOrientation(...
        above, below, lengths, widths, P, image_centerline,...
        linear_convergence_threshold, normalize_homography1D...
    );
end

if display_linear_estimation
    plotProbeReprojection(...
                I, above, below, lengths, widths, P, d, u, X_tip,...
                'Reprojection of linear approximation of probe location'...
            );
end

%% Nonlinear estimation of probe location

[ X_tip, d, u ] = probeTipAndOrientationNonlinear(...
    above, below, lengths, widths, P, X_tip, d,...
    verbose_nonlinear_estimation...
    );

if display_nonlinear_estimation
    plotProbeReprojection(...
                I, above, below, lengths, widths, P, d, u, X_tip,...
                'Reprojection of nonlinearly estimated probe location'...
            );
end