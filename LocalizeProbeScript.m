%% Probe Object 3D Localization
% Use the measured dimensions of the probe object and the locations of
% points along its bands in the image to determine the 3D location of the
% probe object.
%
% ## Usage
%   Modify the parameters, and the paths to input data, in the first code
%   section below, then run.
%
% ## Input
%
% ### Probe detection result
%
% A '.mat' file containing the variables
% 'probe_detection_matches_filtered', and, optionally, 'model_filename',
% output by 'DetectProbeScript.m'. Refer to the documentation of
% 'DetectProbeScript.m' for details.
%
% The probe can be located in 3D space if there are at least three matches
% in 'probe_detection_matches_filtered'. Otherwise, an error will be
% thrown.
%
% ### Probe measurements
%
% The 'model_filename' parameter is the path to a '.mat' file containing a
% structure called 'probe'. 'probe' is a set of user-provided measurements
% of the physical probe. Refer to the documentation of
% 'CreateProbeDetectionModel.m' for details.
%
% This parameter is to be provided separately if not contained in the probe
% detection results described above.
%
% ### Camera projection parameters
%
% The 3 x 4 camera projection matrix, 'P', provided in a '.mat' file.
%
% Note: The camera should be calibrated using the same units of measurement
% as used for the measurements of the probe in the file referred to by
% 'model_filename'.
%
% ### Image containing the probe
%
% The 'I_filename' variable, from the '.mat' file produced by
% 'DetectProbeScript.m', or provided separately, is assumed to contain the
% path to the image in which the probe was detected. This image is used to
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
%   'model_filename' parameter may be provided to this script via the
%   output of 'DetectProbeScript.m' (see above), in which case
%   'model_filename' is actually an indirect parameter of this script. It
%   is copied to the output of this script for completeness.
%
% - 'probe_axis_locations': The `axis_locations` output argument of
%   'localizeProbe()'. Refer to the documentation of 'localizeProbe.m' for
%   details.
%
% - 'probe_axis': The `probe_axis` output argument of 'localizeProbe()'.
%   Refer to the documentation of 'localizeProbe.m' for details.
%
% - 'probe_band_locations': The `band_locations` output argument of
%   'localizeProbe()'. Refer to the documentation of 'localizeProbe.m' for
%   details.
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
model_filename = [];

% Probe detection result
detection_filename = fullfile('..','Data','probeDetectionResults_probePrePaperOcclusion_1_t_rect.mat');

% Camera projection matrix
camera_filename = fullfile('..','Data','20170327_topCameraMatrix_relativeExtrinsics_fixed.mat');

% Image containing the probe
I_filename = [];

% Linear Probe Localization
% Error convergence threshold for linear probe estimation
linear_convergence_threshold = 0.01;
% Normalize lengths when estimating a 1D homography between the probe and
% its image
normalize_homography1D = true;

% Enable or disable nonlinear iterative refinement
enable_nonlinear_estimation = true;

% Debugging tools
verbose.verbose_linear_estimation = false; % Requires `I_filename` to be valid
verbose.display_linear_estimation = false; % Requires `I_filename` to be valid
verbose.verbose_nonlinear_estimation = false;
verbose.display_nonlinear_estimation = false; % Requires `I_filename` to be valid
verbose.display_axis_points = true; % Requires `I_filename` to be valid

%% Load input data
detection_variables_required = {...
        'probe_detection_matches_filtered'...
    };
if ~exist('model_filename', 'var') || isempty(model_filename)
    detection_variables_required = [detection_variables_required 'model_filename'];
end
if ~exist('I_filename', 'var') || isempty(I_filename)
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

%% Pack parameters

params.linear_convergence_threshold = linear_convergence_threshold;
params.normalize_homography1D = normalize_homography1D;
params.enable_nonlinear_estimation = enable_nonlinear_estimation;

%% Location Estimation

[ probe_axis_locations, probe_axis, probe_band_locations ] = localizeProbe(...
    probe, probe_detection_matches_filtered, P, params, verbose, I...
);

%% Save results to a file
save_variables_list = [ parameters_list, {...
        'model_filename',...
        'probe_axis_locations',...
        'probe_axis'...
        'probe_band_locations'...
    } ];
uisave(save_variables_list,'probeLocalizationResults')
