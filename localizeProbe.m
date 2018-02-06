function [ axis_locations, probe_axis, band_locations, hyp ] = localizeProbe(...
       probe, matches_filtered, P, params, varargin...
    )
% LOCALIZEPROBE  Locate the probe object in 3D space
%
% ## Syntax
% axis_locations = localizeProbe(...
%     probe, matches_filtered, P, params [, verbose, I]...
% )
% [ axis_locations, probe_axis ] = localizeProbe(...
%     probe, matches_filtered, P, params [, verbose, I]...
% )
% [ axis_locations, probe_axis, band_locations ] = localizeProbe(...
%     probe, matches_filtered, P, params [, verbose, I]...
% )
% [ axis_locations, probe_axis, band_locations, hyp ] = localizeProbe(...
%     probe, matches_filtered, P, params [, verbose, I]...
% )
%
% ## Description
% axis_locations = localizeProbe(...
%     probe, matches_filtered, P, params [, verbose, I]...
% )
%   Returns the estimated locations of points along the probe's axis of
%   symmetry.
% [ axis_locations, probe_axis ] = localizeProbe(...
%     probe, matches_filtered, P, params [, verbose, I]...
% )
%   Additionally returns the 3D unit direction vector from the tip to the
%   end of the probe.
% [ axis_locations, probe_axis, band_locations ] = localizeProbe(...
%     probe, matches_filtered, P, params [, verbose, I]...
% )
%   Additionally returns a structure describing the reprojection of the
%   probe to image space.
% [ axis_locations, probe_axis, band_locations, hyp ] = localizeProbe(...
%     probe, matches_filtered, P, params [, verbose, I]...
% )
%   Additionally returns the index of the data association hypothesis
%   chosen as a basis for the pose estimate.
%
% ## Input Arguments
%
% probe -- Probe measurements
%   Refer to the documentation of 'CreateProbeDetectionModel.m' for
%   details.
%
% matches_filtered -- Detected and validated edges between coloured bands
%   The `matches_filtered` output argument of 'detectProbe()', describing
%   the points detected on the probe in the image. Refer to the
%   documentation of 'detectProbe.m' for details.
%
% P -- Camera projection matrix
%   The camera matrix (intrinsic and extrinsic parameters). A 3 x 4 matrix.
%
% params -- Fixed parameters
%   Parameters that should be stable across a variety of input images and
%   probe models. `params` is a structure containing the following fields:
%   - linear_convergence_threshold: The convergence threshold used during
%     linear probe location estimation, by the 'homography1D()' function.
%     Refer to the documentation of 'homography1D.m' for details.
%   - normalize_homography1D: Whether or not to normalize point coordinates
%     during the 1D homography estimation subroutine of linear probe
%     location estimation, 'homography1D()'. Refer to the documentation of
%     'homography1D.m' for details.
%   - enable_nonlinear_estimation: If `true`, follow linear probe location
%     estimation with an iterative optimization procedure that refines the
%     location estimate.
%
% verbose -- Debugging flags
%   If recognized fields of `verbose` are true, corresponding graphical or
%   console output will be generated for debugging purposes. `I` must be
%   passed if any debugging flags that trigger graphical output are set.
%
%   All debugging flags default to false if `verbose` is not passed.
%
% I -- Image containing the probe
%   An RGB image to use as a background for an overlay of the probe in the
%   context of graphical debugging output.
%
% ## Output Arguments
%
% axis_locations -- Points along the probe axis
%   A structure vector describing the estimated positions of the points
%   along the central axis of the probe object. `band_locations` is useful
%   for assessing reprojection error, whereas `axis_locations` provides
%   data for further applications (e.g. producing a point cloud by
%   localizing the probe tip in a sequence of probe positions). Each
%   element contains the locations of a point on the line through the
%   centre of the probe, corresponding to a length measurement in
%   `probe.lengths`. The points are in order, so the first point is the tip
%   of the probe, while the last point is the back end of the probe.
%
%   The fields of the structure vector are as follows:
%   - imagePoint: A two-element row vector containing the estimated
%       position of the point in the image.
%   - objectPoint: A three-element row vector containing the estimated
%       position of the point in 3D space. (The point is in the coordinate
%       frame determined by the extrinsic parameters of the input camera
%       matrix, `P`.)
%
% probe_axis -- Probe direction vector
%   A 3D unit row vector containing the estimated direction of the probe
%   object in space. 'probe_axis' points from the tip of the probe
%   (corresponding to the first measurement in `probe.lengths`) to the end
%   of the probe.
%
% band_locations -- Points along the edges of the probe
%   A structure vector describing the positions of the detected interest
%   points on the probe. Each element contains the locations of the two
%   endpoints of one detected edge between probe bands. This structure is a
%   version of a cell of `matches_filtered`, `matches_filtered(hyp)`,
%   wherein points are reprojected from their estimated 3D locations,
%   instead of detected in the image. In contrast to `axis_locations`,
%   which contains data for all edges between probe bands, `band_locations`
%   only contains data for the edges which were detected in the image.
%
%   The fields of the structure vector are as follows:
%   - index: The index of the detected edge between two bands of colour
%       on the probe. The indices of edges correspond to the order of the
%       edges along the major axis of the probe in the image, from left to
%       right.
%   - pointAbovePCAMajorAxis: A two-element row vector containing the
%       reprojected position of the point where the edge meets the lower
%       border of the probe in the image. (This is the point "above" the
%       axis of the probe in the sense that it has a pixel y-coordinate
%       that is larger than that of the axis at the same pixel
%       x-coordinate.)
%   - pointBelowPCAMajorAxis: Similar to 'pointAbovePCAMajorAxis', but
%       contains the reprojected position of the point where the edge meets
%       the upper border of the probe in the image.
%   - matchedLengthIndex: The index of the edge between coloured bands of
%       the physical probe that is matched with the detected edge.
%       Specifically, this is the index into the user-provided measurements
%       of the probe, `probe.lengths` and `probe.widths`.
%
% hyp -- Data association hypothesis index
%   An integer specifying which cell in `matches_filtered` was selected as
%   the data association hypothesis used for pose estimation.
%
% ## Algorithm
%
% Probe location is a two-step process: Linear initialization, followed by
% nonlinear refinement. The initial estimate of the probe's location is
% obtained by first approximating the axis of symmetry of the probe with
% the centerline passing through each of the detected probe band edges.
% Then, 'probeTipAndOrientation()' solves a linear system to find the
% location of the probe which minimizes algebraic error. In the second
% step, 'probeTipAndOrientationNonlinear()' refines the initial solution by
% minimizing reprojection error using the Levenberg-Marquardt algorithm.
%
% The linear initialization is performed for each data association
% hypothesis in `matches_filtered` (i.e. for each cell in
% `matches_filtered`). The hypothesis which results in the lowest
% reprojection error is selected for the nonlinear refinement step, and the
% index of this chosen hypothesis is output in `hyp`.
%
% ## Notes
% - The camera calibration matrix should be expressed using the same units
%   of measurement as the values in `probe`.
% - The probe can be located in 3D space if there are at least three
%   elements (i.e. correspondences between measured and detected probe band
%   edges) in `matches_filtered`. Otherwise, an error will be thrown.
% - Presently, all detected points are assumed to be valid, and correctly
%   associated with measured points on the probe.
%
% See also detectProbe, pca, probeTipAndOrientation, probeTipAndOrientationNonlinear

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created January 10, 2018

%% Parse input arguments

nargoutchk(1,4);
narginchk(4,6);

if ~isempty(varargin)
    verbose = varargin{1};
    verbose_linear_estimation = verbose.verbose_linear_estimation;
    display_linear_estimation = verbose.display_linear_estimation;
    verbose_nonlinear_estimation = verbose.verbose_nonlinear_estimation;
    display_nonlinear_estimation = verbose.display_nonlinear_estimation;
    display_axis_points = verbose.display_axis_points;
else
    verbose_linear_estimation = false;
    display_linear_estimation = false;
    verbose_nonlinear_estimation = false;
    display_nonlinear_estimation = false;
    display_axis_points = false; 
end

if verbose_linear_estimation || display_linear_estimation ||...
   display_nonlinear_estimation || display_axis_points
    if length(varargin) < 2
        error('An image must be passed to produce the requested debugging visualizations.');
    else
        I = varargin{2};
    end
end

% Only keep sufficiently long series of matches
n_hypotheses = length(matches_filtered);
hypothesis_filter = false(n_hypotheses, 1);
for hyp = 1:n_hypotheses
    hypothesis_filter(hyp) = (length(matches_filtered{hyp}) >= 3);
end
if ~any(hypothesis_filter)
    error('Insufficient number of probe band edge correspondences for probe location estimation.');
end

%% Linear estimation of probe location

% Pick the best hypothesis, in terms of sum of squared reprojection error
min_error = Inf;
hyp_final = 0;
above = cell(n_hypotheses, 1);
below = cell(n_hypotheses, 1);
lengths = cell(n_hypotheses, 1);
widths = cell(n_hypotheses, 1);
X_tip_per_hypothesis = zeros(n_hypotheses, 3);
probe_axis_per_hypothesis = zeros(n_hypotheses, 3);
for hyp = 1:n_hypotheses
    if hypothesis_filter(hyp)
        hypothesis = matches_filtered{hyp};
        above{hyp} = vertcat(hypothesis(:).pointAbovePCAMajorAxis);
        below{hyp} = vertcat(hypothesis(:).pointBelowPCAMajorAxis);
        lengths{hyp} = vertcat(hypothesis(:).matchedLength);
        widths{hyp} = vertcat(hypothesis(:).matchedWidth);

        % Estimated centerline of the probe in the image
        allPoints = [ above{hyp}; below{hyp} ];
        [ coeff, ~, ~, ~, ~, mu ] = pca(allPoints);

        % Express the PCA component vectors as lines in the space of the original data
        axes = pcaAxes2D( coeff, mu );

        % The first axis is the estimated centerline of the probe
        %
        % I have not found a linear method for estimating the true axis of the
        % probe (taking projective distortion into account) that remains
        % numerically stable as the width of the probe decreases. However, as the
        % width of the probe decreases, the first PCA axis becomes an increasingly
        % good approximation of the true probe axis.
        image_centerline = axes(1, :);

        if verbose_linear_estimation
            [...
                X_tip_per_hypothesis(hyp, :), probe_axis_per_hypothesis(hyp, :)...
            ] = probeTipAndOrientation(...
                above{hyp}, below{hyp}, lengths{hyp}, widths{hyp}, P, image_centerline,...
                params.linear_convergence_threshold, params.normalize_homography1D, I...
            );
        else
            [...
                X_tip_per_hypothesis(hyp, :), probe_axis_per_hypothesis(hyp, :)...
            ] = probeTipAndOrientation(...
                above{hyp}, below{hyp}, lengths{hyp}, widths{hyp}, P, image_centerline,...
                params.linear_convergence_threshold, params.normalize_homography1D...
            );
        end
        
        % Check error
        [above_reprojected, below_reprojected] = reprojectProbe(...
            lengths{hyp}, widths{hyp}, P, probe_axis_per_hypothesis(hyp, :), X_tip_per_hypothesis(hyp, :)...
        );
        reprojected = [above_reprojected; below_reprojected];
        err = allPoints - reprojected;
        err = sum(sum(err.^2)) / size(allPoints, 1);
        if err < min_error
            min_error = err;
            hyp_final = hyp;
        end
    end
end

hyp = hyp_final;
above = above{hyp};
below = below{hyp};
lengths = lengths{hyp};
widths = widths{hyp};
probe_axis = probe_axis_per_hypothesis(hyp, :);
X_tip = X_tip_per_hypothesis(hyp, :);

if display_linear_estimation
    plotProbeReprojection(...
                I, [above; below], lengths, widths, P, probe_axis, X_tip,...
                'Reprojection of linear approximation of probe location'...
            );
end

%% Nonlinear estimation of probe location

if params.enable_nonlinear_estimation
    [ X_tip, probe_axis ] = probeTipAndOrientationNonlinear(...
        above, below, lengths, widths, P, X_tip, probe_axis,...
        verbose_nonlinear_estimation...
        );

    if display_nonlinear_estimation
        plotProbeReprojection(...
                    I, [above; below], lengths, widths, P, probe_axis, X_tip,...
                    'Reprojection of nonlinearly estimated probe location'...
                );
    end
end

%% Output

n = size(probe.lengths, 1);
axis_points = repmat(X_tip, n, 1) + repmat(probe.lengths, 1, 3) .* repmat(probe_axis, n, 1);
axis_points_image = (P * [axis_points, ones(n,1)].').';
axis_points_image = axis_points_image(:, 1:2) ./ repmat(axis_points_image(:, end), 1, 2);

if display_axis_points
    figure;
    imshow(I);
    hold on
    scatter(axis_points_image(:, 1), axis_points_image(:, 2), 'ro');
    hold off
    title('Points along the probe axis')
end

axis_locations = struct(...
        'imagePoint', num2cell(axis_points_image, 2),...
        'objectPoint', num2cell(axis_points, 2)...
    );

if nargout > 2
    [above_reprojected, below_reprojected] = reprojectProbe(...
        lengths, widths, P, probe_axis, X_tip...
        );

    band_locations = struct(...
            'index', {matches_filtered{hyp}.index}.',...
            'pointAbovePCAMajorAxis', num2cell(above_reprojected, 2),...
            'pointBelowPCAMajorAxis', num2cell(below_reprojected, 2),...
            'matchedLengthIndex', {matches_filtered{hyp}.matchedLengthIndex}.'...
        );
end

end
