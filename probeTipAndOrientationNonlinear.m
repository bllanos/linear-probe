function [ X_tip, d, u ] = probeTipAndOrientationNonlinear( above, below, lengths, widths, P, X_tip, d, varargin )
%PROBETIPANDORIENTATIONNONLINEAR Nonlinear estimation of probe tip and orientation
%
% ## Syntax
% [X_tip, d, u] = probeTipAndOrientationNonlinear(...
%     above, below, lengths, widths, P, X_tip, d [, verbose]
% )
%
% ## Description
% [X_tip, d, u] = probeTipAndOrientationNonlinear(...
%     above, below, lengths, widths, P, X_tip, d [, verbose]
% )
%   Uses a nonlinear iterative method (minimization of geometric error)
%   to calculate the position of the probe tip and the orientation of the
%   probe, given the positions of points on the probe in the image, and
%   given calibration information for the camera, and the probe dimensions.
%
% ## Input Arguments
%
% above -- Interest points along lower edge of probe
%   Points located at higher y-coordinates than the probe midline in the
%   image. An n x 2 array of image coordinates.
%
% below -- Interest points along upper edge of probe
%   Points located at lower y-coordinates than the probe midline in the
%   image. An n x 2 array of image coordinates. `below(i, :)` is the point
%   oppposite `above(i, :)` on an edge between two colour bands of the
%   probe.
%
% lengths -- Probe length measurements
%   A vector of length `n` containing the measured physical distances of
%   the pairs of points in `above` and `below` from the tip of the probe.
%
% widths -- Probe diameter measurements
%   A vector of length `n` containing the measured diameters of the probe
%   at each of the locations corresponding to the points in `above` and
%   `below`.
%
% P -- Camera projection matrix
%   The camera matrix (intrinsic and extrinsic parameters) corresponding to
%   the image points in `above` and `below`. The camera must be a finite
%   camera, not an affine camera.
%
% X_tip -- Probe tip guess
%   An initial guess for the 3D position of the probe tip. A 3 x 1 vector.
%
% d -- Probe orientation vector guess
%   An initial guess for the 1 x 3 unit vector aligned with the axis of the
%   probe in space.
%
% verbose -- Debugging flag
%   If true, console output will be generated for debugging purposes.
%
%   Defaults to false if not passed.
%
% ## Output Arguments
%
% X_tip -- Probe tip
%   The 3D position of the probe tip, output as a 1 x 3 vector.
%
% d -- Probe orientation vector
%   A 1 x 3 unit vector aligned with the axis of the probe in space.
%
% u -- Probe normal vector
%   A 1 x 3 unit vector perpendicular to the plane formed by the axis of
%   the probe in space and the camera centre, pointing towards the bottom
%   of the image.
%
% ## References
% - R. Hartley and A. Zisserman. Multiple View Geometry in Computer Vision,
%   2nd Edition. Cambridge, UK: Cambridge University Press, 2003.
%
% See also probeTipAndOrientation, lsqnonlin

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created March 8, 2017

nargoutchk(1, 3);
narginchk(7, 8);

if ~isempty(varargin)
    verbose = varargin{1};
else
    verbose = false;
end

% Initial Guess
X_end = X_tip + d * max(lengths);
p0 = zeros(1, 6);
if verbose
    disp('Initial guess:')
    X_tip %#ok<NOPRT>
    X_end %#ok<NOPRT>
    disp('Perturbation (the actual initial guess for the optimizer):')
    disp(p0)
end

allPoints = [above; below];

    function [X_tip, X_end] = parseSolution( p, X_tip_0, X_end_0 )
        % p = [ X_tip perturbation, X_end perturbation ]
        X_tip = X_tip_0 + p(1:3);
        X_end = X_end_0 + p(4:end);
    end

    function [ error ] = costFunction( p )
        [X_tip_i, X_end_i] = parseSolution( p, X_tip, X_end );
        d_i = probeOrientationFromEndpoints( P, X_tip_i, X_end_i );
        [above_reprojected, below_reprojected] = reprojectProbe(...
            lengths, widths, P, d_i, X_tip_i...
        );
        reprojected = [above_reprojected; below_reprojected];
        error = allPoints - reprojected;
        error = sqrt(sum(error.^2, 2));
        % I could compute the Jacobian to speed up the optimization process
    end

% Invoke optimizer
options = optimoptions('lsqnonlin',...
    'FunValCheck', 'on', 'Algorithm', 'levenberg-marquardt'...
);
if verbose
    options.Display = 'iter-detailed';
else
    options.Display = 'none';
end
p = lsqnonlin(@costFunction, p0, [], [], options);

[X_tip, X_end] = parseSolution( p, X_tip, X_end );
[d, u] = probeOrientationFromEndpoints( P, X_tip, X_end );

if verbose
    disp('Final guess from optimizer:')
    disp(p)
    X_tip %#ok<NOPRT>
    X_end %#ok<NOPRT>
    d %#ok<NOPRT>
    u %#ok<NOPRT>
end

end

