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
% quantified by calculating image reprojection error.
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
%   - gt: A four-element row vector. The first two elements contain the
%       user-marked position of the interest point on the lower border of
%       the probe in the image. The last two elements contain the
%       user-marked position of the interest point on the upper border of
%       the probe in the image.
%   - reprojected: A four-element row vector containing the reprojected
%       positions of the interest points on the lower (first two elements)
%       and upper (last two elements) border of the probe in the image.
%   - error: A two-element row vector containing The Euclidean
%       distances between the reprojected and ground-truth points on the
%       lower, and upper, borders of the probe in the image, respectively.
%
%   Ground truth information for the tapered tip(s) of the probe is
%   duplicated so that both the 'gt' and 'error' fields have values that
%   are the same sizes as those for the other points. The 'index' field can
%   be used to distinguish probe tips from other contour points.
%
% - 'error_mean_2d': The average of the 'error' field values in
%   'visible_points_error' (avoiding double-counting the probe tips).
%
% - 'error_max_2d': The maximum of the 'error' field values in
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

%% Input data and parameters

% List of parameters to save with results
parameters_list = {
        'localization_filename',...
        'visible_points_filename'...
        'point_alignment_outlier_threshold'
    };
% The following variables are loaded from the file pointed to by
% 'localization_filename', unless they are defined below:
parameters_list = [parameters_list {
        'model_filename',...
        'camera_filename',...
        'I_filename'
    }];

% Probe localization results
localization_filename = 'C:\Users\llanos\Google Drive\PointProbing\Data and results\20160811_bambooSkewerProbe\20170409_bambooSkewer_orangeBlue_probeLocalizationResults_topCamera_rect.mat';

% Measurements of the probe object
model_filename = [];

% Camera projection matrix
camera_filename =  'C:\Users\llanos\Google Drive\PointProbing\Data and results\20160811_bambooSkewerProbe\camera_calibration_from_20160609\20170305_bottomCameraMatrix_identityExtrinsics.mat';

% Image containing the probe
I_filename = 'C:\Users\llanos\Google Drive\PointProbing\Data and results\20160811_bambooSkewerProbe\undistorted\probePrePaperOcclusion_1_b_rect.bmp';

% Ground truth points
% If none is provided, set it to an empty array (`[]`).
visible_points_filename = 'C:\Users\llanos\Google Drive\PointProbing\Data and results\20160811_bambooSkewerProbe\ground_truth\probePrePaperOcclusion_1_b_groundTruth.mat';

% Parameters for interpreting annotated points
point_alignment_outlier_threshold = 5;

% Flags to control verbose or graphical output
display_reprojection = true;
display_model_from_visible_points = true;
display_reprojection_statistics = true;

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

%% Assess reprojection error

if single_view_gt_available
    
    % Classify marked points
    [ ~, ~, ~, visible_points_model ] = bilateralModel(...
        visible_points(:, 2:3), point_alignment_outlier_threshold, true...
    );
    if display_model_from_visible_points
        fg = figure; %#ok<UNRCH>
        imshow(I);
        plotBilateralModel( visible_points_model, [], [], [], fg);
        title('Classified annotated points (blue, black = tips; green, red = above/below first PCA component; yellow = unmatched)');
    end
    if ~isempty(visible_points_model.unmatched)
        error(sprintf(['Not all probe colour band junctions are marked with two points - One on each edge of the probe.\n',...
            'Consider increasing the outlier detection threshold used when pairing points.'])); %#ok<SPERR>
    end
    
    % Associate with reprojected points
    n_visible_points = size(visible_points, 1);
    n_visible_points_pairs = size(visible_points_model.above, 1);
    visible_length_indices = zeros(n_visible_points_pairs, 1);
    for i = 1:n_visible_points_pairs
        visible_length_indices(i) =...
            visible_points(...
                all(...
                    repmat(visible_points_model.above(i, :), n_visible_points, 1)...
                    == visible_points(:, 2:end)...
                    , 2 ...
                ),...
                1 ...
            );
    end
    
    % Make sure the points are sorted by band index (in order for the
    % output to be better organized)
    [visible_length_indices, sorting_map] = sort(visible_length_indices);
    above_visible = visible_points_model.above(sorting_map, :);
    below_visible = visible_points_model.below(sorting_map, :);
    
    % Obtain reprojected points for the interior of the probe
    [above_reprojected, below_reprojected] = reprojectProbe(...
        probe.lengths(visible_length_indices),...
        probe.widths(visible_length_indices), P, d, X_tip...
    );

    % Calculate reprojection error
    gt = [above_visible, below_visible];
    reprojected = [above_reprojected, below_reprojected];
    error_interior = gt - reprojected;
    error_interior = error_interior.^2;
    error_interior = [ sqrt(sum(error_interior(:, 1:2), 2)), sqrt(sum(error_interior(:, 3:end), 2)) ];
    
    visible_points_error = struct(...
        'index', num2cell(visible_length_indices),...
        'gt', num2cell(gt, 2),...
        'reprojected', num2cell(reprojected, 2),...
        'error', num2cell(error_interior, 2)...
    );

    error_sum = sum(sum(error_interior));
    error_max_2d = max(max(error_interior));
    error_count = numel(error_interior);
    
    % Repeat for the probe tips
    for name = {'head', 'tail'}
        if isfield(visible_points_model, name{1})
            tip_index = visible_points(...
                    all(...
                        repmat(visible_points_model.(name{1}), n_visible_points, 1)...
                        == visible_points(:, 2:end)...
                        , 2 ...
                    ),...
                    1 ...
                );
            tip_reprojected = reprojectProbe(...
                probe.lengths(tip_index),...
                probe.widths(tip_index), P, d, X_tip...
            );
            tip_error = sqrt(sum((visible_points_model.(name{1}) - tip_reprojected).^2));
            tip_struct = struct(...
                'index', tip_index,...
                'gt', repmat(visible_points_model.(name{1}), 1, 2),...
                'reprojected', repmat(tip_reprojected, 1, 2),...
                'error', repmat(tip_error, 1, 2)...
            );
            error_sum = error_sum + tip_error;
            error_max_2d = max(error_max_2d, tip_error);
            error_count = error_count + 1;

            if tip_index < visible_length_indices(1)
                visible_points_error = [ tip_struct; visible_points_error ]; %#ok<AGROW>
            else
                visible_points_error = [ visible_points_error; tip_struct ]; %#ok<AGROW>
            end
        end
    end
    
    error_mean_2d = error_sum / error_count;
    
    if display_reprojection_statistics
        disp('Reprojection error:') %#ok<UNRCH>
        visible_points_error_display = struct2table(visible_points_error);
        disp(visible_points_error_display);
        fprintf('Mean reprojection error: %g\n', error_mean_2d);
        fprintf('Maximum reprojection error: %g\n', error_max_2d);
        
        figure;
        error_mean_pairwise = mean(vertcat(visible_points_error(:).error), 2);
        plot(error_mean_pairwise);
        title('Reprojection error averaged over points above and below the probe axis')
        xlabel('Index of location along probe')
        ylabel('Reprojection error [pixels]')
    end
end


%% Save results to a file
if single_view_gt_available
    save_variables_list = [ parameters_list, {...
            'visible_points_error',...
            'error_mean_2d',...
            'error_max_2d'
        } ];
    uisave(save_variables_list,'probeLocalizationEvaluation')
end