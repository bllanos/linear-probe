function [X_tip, d, u] = probeTipAndOrientation( above, below, lengths, widths, P, image_line, threshold, normalize, varargin )
%PROBETIPANDORIENTATION Linear estimation of probe tip and orientation
%
% ## Syntax
% [X_tip, d, u] = probeTipAndOrientation(...
%     above, below, lengths, widths, P, image_line, threshold, normalize...
%     [, I]
% )
%
% ## Description
% [X_tip, d, u] = probeTipAndOrientation(...
%     above, below, lengths, widths, P, image_line, threshold, normalize...
%     [, I]
% )
%   Uses a linear method (minimization of algebraic, not geometric error)
%   to calculate the position of the probe tip and the orientation of the
%   probe, given the positions of points on the probe in the image, and
%   given calibration information for the camera and the probe dimensions.
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
% image_line -- Estimated probe axis in image
%   A 3-vector containing the homogenous coordinates of a line in the image
%   approximating the probe axis.
%
% threshold -- Convergence threshold
%   When the relative change in the residual error between iterations of
%   the linear solution drops below this threshold, return the current
%   solution. (Also used as the convergence threshold for the
%   'homography1D.m' subroutine.)
%
% normalize -- Choice of image point normalization
%   Whether or not to perform normalization of point coordinates
%   prior to estimating a 1D homography as a subroutine, 'homography1D.m'.
%
% I -- Image for debugging visualization
%   If passed, the function will plot intermediate results on top of this
%   image, for debugging/visualization purposes. It will also output
%   debugging information to the console.
%
% ## Output Arguments
%
% X_tip -- Probe tip
%   The 3D position of the probe tip, output as a 3 x 1 vector.
%
% d -- Probe orientation vector
%   A 3 x 1 unit vector aligned with the axis of the probe in space.
%
% u -- Probe normal vector
%   A 3 x 1 unit vector perpendicular to the plane formed by the axis of
%   the probe in space and the camera centre, pointing towards the bottom
%   of the image.
%
% ## Side Effects
%
% - Opens figures for debugging visualizations, and produces console
%   output, if `I` is passed as an input argument.
%
% ## References
% - R. Hartley and A. Zisserman. Multiple View Geometry in Computer Vision,
%   2nd Edition. Cambridge, UK: Cambridge University Press, 2003.
%
% See also homography1D

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created March 8, 2017

nargoutchk(3, 3);
narginchk(8, 9);

verbose = ~isempty(varargin);
if verbose
    I = varargin{1};
    image_size = size(I);
    image_size = image_size(1:2);
end
i = 0; % Iteration counter

% Find an initial estimate of `u`
    function [u] = uFromImageLine(image_line)
        % The plane back-projected through the centerline of the probe
        % Result 8.2 from Multiple View Geometry in Computer Vision, 2nd
        % Edition
        centerline_plane = image_line * P;
        
        % Normal vector to the plane
        u = centerline_plane(1:3);
        u = u ./ repmat(norm(u), 1, 3); % Normalize
        % Make sure the vector points to the bottom of the image
        u_image = (P * [u 0].').';
        u_image = u_image(1:2) ./ repmat(u_image(end), 1, 2);
        y_centerline = -(image_line(1) * u_image(1) + image_line(3)) / image_line(2);
        if y_centerline > u_image(2)
            u = -u;
        end
    end

u = uFromImageLine(image_line);

n = size(above, 1);
nAll = 2 * n;
above = [above, ones(n, 1)];
below = [below, ones(n, 1)];
allPoints = [above; below];

l = repmat(lengths, 2, 1);

% Find an initial estimate of the probe tip in the image
[camera_xAxis, camera_yAxis, camera_zAxis] = cameraAxes( P );

    % `X_tip_image` is expressed in normalized homogenous coordinates
    function [X_tip_image, X_end_image, tangent_3D] = parametersFromMidline(image_line)
        % Convert the points to 1D coordinates on the estimated line
        image_line_points = closestPointOnLine(repmat(image_line, nAll, 1), allPoints(:, 1:2));
        tangent = [image_line(2), -image_line(1)];
        tangent = tangent / norm(tangent);
        image_distances = image_line_points - repmat(image_line_points(1, :), nAll, 1);
        image_distances = dot(image_distances, repmat(tangent, nAll, 1), 2);
        H = homography1D( l, image_distances, normalize, threshold );
        X_tip_image = (H * [0; 1]).';
        X_tip_image = X_tip_image(1) / X_tip_image(2);
        X_tip_image = [image_line_points(1, :) + tangent * X_tip_image, 1];
        X_end_image = (H * [max(lengths); 1]).';
        X_end_image = X_end_image(1) / X_end_image(2);
        X_end_image = [image_line_points(1, :) + tangent * X_end_image, 1];
        
        if verbose
            figure;
            imshow(I);
            hold on
            line_points_plotting = lineToBorderPoints(image_line, image_size);
            line(line_points_plotting([1,3]), line_points_plotting([2,4]), 'Color', 'c');
            scatter(allPoints(:, 1), allPoints(:, 2), 'g.');
            scatter(image_line_points(:, 1), image_line_points(:, 2), 'bo');
            reprojected_points = (H * [l ones(nAll, 1)].').';
            reprojected_points = reprojected_points(:, 1) ./ reprojected_points(:, 2);
            reprojected_points = image_line_points(1, :) + repmat(tangent, nAll, 1) .* repmat(reprojected_points, 1, 2);
            scatter(reprojected_points(:, 1), reprojected_points(:, 2), 'r.');
            scatter(X_tip_image(1), X_tip_image(2), 'm+');
            scatter(X_end_image(1), X_end_image(2), 'c+');
            hold off
            legend(...
                'Estimated axis', 'Detected points', 'Projected onto axis',...
                'Reprojected from homography', 'Estimated tip',...
                'Estimated furthest detected point'...
                );
            title(sprintf('1D Homography Estimation for iteration %d', i))
        end
        
        % Find the tangent vector in world coordinates (aligned with the
        % image plane)
        tangent_3D = camera_xAxis * tangent(1) + camera_yAxis * tangent(2);
        tangent_3D = tangent_3D(1:3, :);
    end

[X_tip_image, X_end_image, tangent_3D] = parametersFromMidline(image_line);

i = 1;

% Parameterize the probe tip as a point on the ray through `X_tip_image`,
% with some offset in the direction of `tangent_3D`.
P_inv = pinv(P); % Pseudoinverse
P_center = null(P).'; % Camera center
% Assume P_center(end) ~= 0 (i.e. Finite camera)
P_center = P_center ./ repmat(P_center(end), 1, 4);
X_tip_basis_ray = (P_inv * X_tip_image.').';
X_end_basis_ray = (P_inv * X_end_image.').';

% From P * (X_tip + l_i * [d; 0] + r * [u; 0]) ~ x_i
% So cross(P * (X_tip + l_i * [d; 0] + r * [u; 0]), x_i) = 0
% An image line and its known cross ratios have 5 degrees of freedom, so it
% is necessary to reduce the number of parameters:
%   Substitute `X_tip = P_center + lambda1 * X_tip_basis_ray + lambda2 * [tangent_3D; 0]`
%   Note: `X_tip = P_center + lambda1 * X_tip_basis_ray + lambda2 * [u; 0]`
%     produces a sub-optimal solution, in which the probe axis tilts until
%     it is a diagonal of the approximate rectangle formed by the points
%     detected along the probe.
%
% Normalization of `d` is a nonlinear constraint.
%   Reparameterize as P * (X_tip * (1 - k) + X_end * k + r * [u; 0]) ~ x_i,
%   Where k is the normalized version of 1_i (where the last point
%   detected on the probe has k = 1).
%   `X_end = P_center + lambda3 * X_end_basis_ray + lambda4 * [tangent_3D; 0]`
%
% Express as a matrix A * [lambda1; lambda2; lambda3; lambda4] = b
x1 = allPoints(:, 1);
x2 = allPoints(:, 2);
r = repmat(widths / 2, 2, 1); % Take radii, not diameters
r(n+1:end) = -r(n+1:end); % Account for the opposition between `above` and `below`

k = l / max(lengths);

P1_1 = P(1,1);
P1_2 = P(1,2);
P1_3 = P(1,3);
P1_4 = P(1,4);
P2_1 = P(2,1);
P2_2 = P(2,2);
P2_3 = P(2,3);
P2_4 = P(2,4);
P3_1 = P(3,1);
P3_2 = P(3,2);
P3_3 = P(3,3);
P3_4 = P(3,4);
P_center1 = P_center(1);
P_center2 = P_center(2);
P_center3 = P_center(3);
P_center4 = P_center(4);

    function A = rhs(X_tip_basis_ray, X_end_basis_ray, tangent_3D)
        Xb1 = X_tip_basis_ray(1);
        Xb2 = X_tip_basis_ray(2);
        Xb3 = X_tip_basis_ray(3);
        Xb4 = X_tip_basis_ray(4);
        Xe1 = X_end_basis_ray(1);
        Xe2 = X_end_basis_ray(2);
        Xe3 = X_end_basis_ray(3);
        Xe4 = X_end_basis_ray(4);
        t1 = tangent_3D(1);
        t2 = tangent_3D(2);
        t3 = tangent_3D(3);
        A = [
            (P2_1.*Xb1.*(k - 1) - x2.*(P3_1.*Xb1.*(k - 1) + P3_2.*Xb2.*(k - 1) + P3_3.*Xb3.*(k - 1) + P3_4.*Xb4.*(k - 1)) + P2_2.*Xb2.*(k - 1) + P2_3.*Xb3.*(k - 1) + P2_4.*Xb4.*(k - 1)),...
            (P2_1.*t1.*(k - 1) - x2.*(P3_1.*t1.*(k - 1) + P3_2.*t2.*(k - 1) + P3_3.*t3.*(k - 1)) + P2_2.*t2.*(k - 1) + P2_3.*t3.*(k - 1)),...
            (x2.*(P3_1.*Xe1.*k + P3_2.*Xe2.*k + P3_3.*Xe3.*k + P3_4.*Xe4.*k) - P2_1.*Xe1.*k - P2_2.*Xe2.*k - P2_3.*Xe3.*k - P2_4.*Xe4.*k),...
            (x2.*(P3_1.*k.*t1 + P3_2.*k.*t2 + P3_3.*k.*t3) - P2_1.*k.*t1 - P2_2.*k.*t2 - P2_3.*k.*t3);
            
            (x1.*(P3_1.*Xb1.*(k - 1) + P3_2.*Xb2.*(k - 1) + P3_3.*Xb3.*(k - 1) + P3_4.*Xb4.*(k - 1)) - P1_1.*Xb1.*(k - 1) - P1_2.*Xb2.*(k - 1) - P1_3.*Xb3.*(k - 1) - P1_4.*Xb4.*(k - 1)),...
            (x1.*(P3_1.*t1.*(k - 1) + P3_2.*t2.*(k - 1) + P3_3.*t3.*(k - 1)) - P1_1.*t1.*(k - 1) - P1_2.*t2.*(k - 1) - P1_3.*t3.*(k - 1)),...
            (P1_1.*Xe1.*k - x1.*(P3_1.*Xe1.*k + P3_2.*Xe2.*k + P3_3.*Xe3.*k + P3_4.*Xe4.*k) + P1_2.*Xe2.*k + P1_3.*Xe3.*k + P1_4.*Xe4.*k),...
            (P1_1.*k.*t1 - x1.*(P3_1.*k.*t1 + P3_2.*k.*t2 + P3_3.*k.*t3) + P1_2.*k.*t2 + P1_3.*k.*t3);
            
            (x2.*(P1_1.*Xb1.*(k - 1) + P1_2.*Xb2.*(k - 1) + P1_3.*Xb3.*(k - 1) + P1_4.*Xb4.*(k - 1)) - x1.*(P2_1.*Xb1.*(k - 1) + P2_2.*Xb2.*(k - 1) + P2_3.*Xb3.*(k - 1) + P2_4.*Xb4.*(k - 1))),...
            (x2.*(P1_1.*t1.*(k - 1) + P1_2.*t2.*(k - 1) + P1_3.*t3.*(k - 1)) - x1.*(P2_1.*t1.*(k - 1) + P2_2.*t2.*(k - 1) + P2_3.*t3.*(k - 1))),...
            (x1.*(P2_1.*Xe1.*k + P2_2.*Xe2.*k + P2_3.*Xe3.*k + P2_4.*Xe4.*k) - x2.*(P1_1.*Xe1.*k + P1_2.*Xe2.*k + P1_3.*Xe3.*k + P1_4.*Xe4.*k)),...
            (x1.*(P2_1.*k.*t1 + P2_2.*k.*t2 + P2_3.*k.*t3) - x2.*(P1_1.*k.*t1 + P1_2.*k.*t2 + P1_3.*k.*t3))
            ];
    end

A = rhs(X_tip_basis_ray, X_end_basis_ray, tangent_3D);

    function b = lhs(u)
        u1 = u(1);
        u2 = u(2);
        u3 = u(3);
        b = [
            x2.*(P3_4.*(P_center4.*k - P_center4.*(k - 1)) + P3_1.*(P_center1.*k + r.*u1 - P_center1.*(k - 1)) + P3_2.*(P_center2.*k + r.*u2 - P_center2.*(k - 1)) + P3_3.*(P_center3.*k + r.*u3 - P_center3.*(k - 1))) - P2_4.*(P_center4.*k - P_center4.*(k - 1)) - P2_1.*(P_center1.*k + r.*u1 - P_center1.*(k - 1)) - P2_2.*(P_center2.*k + r.*u2 - P_center2.*(k - 1)) - P2_3.*(P_center3.*k + r.*u3 - P_center3.*(k - 1));
            P1_4.*(P_center4.*k - P_center4.*(k - 1)) - x1.*(P3_4.*(P_center4.*k - P_center4.*(k - 1)) + P3_1.*(P_center1.*k + r.*u1 - P_center1.*(k - 1)) + P3_2.*(P_center2.*k + r.*u2 - P_center2.*(k - 1)) + P3_3.*(P_center3.*k + r.*u3 - P_center3.*(k - 1))) + P1_1.*(P_center1.*k + r.*u1 - P_center1.*(k - 1)) + P1_2.*(P_center2.*k + r.*u2 - P_center2.*(k - 1)) + P1_3.*(P_center3.*k + r.*u3 - P_center3.*(k - 1));
            x1.*(P2_4.*(P_center4.*k - P_center4.*(k - 1)) + P2_1.*(P_center1.*k + r.*u1 - P_center1.*(k - 1)) + P2_2.*(P_center2.*k + r.*u2 - P_center2.*(k - 1)) + P2_3.*(P_center3.*k + r.*u3 - P_center3.*(k - 1))) - x2.*(P1_4.*(P_center4.*k - P_center4.*(k - 1)) + P1_1.*(P_center1.*k + r.*u1 - P_center1.*(k - 1)) + P1_2.*(P_center2.*k + r.*u2 - P_center2.*(k - 1)) + P1_3.*(P_center3.*k + r.*u3 - P_center3.*(k - 1)))
       ];
    end

b = lhs(u);

% Minimize L2 norm of (A.p - b)
p = A \ b;

% Update the estimate of `u`
    function updateSolution(p)
        X_tip = P_center + p(1) * X_tip_basis_ray + p(2) * [tangent_3D.' 0];
        X_end = P_center + p(3) * X_end_basis_ray + p(4) * [tangent_3D.' 0];
        
        d = X_end(1:3) - X_tip(1:3);
        estimated_length = norm(d);
        scale = max(lengths) / estimated_length;
        
        % Rescale so that the probe has the correct length
        X_tip = P_center + scale * (p(1) * X_tip_basis_ray + p(2) * [tangent_3D.' 0]);
        X_end = P_center + scale * (p(3) * X_end_basis_ray + p(4) * [tangent_3D.' 0]);
        
        X_tip_image = (P * X_tip.').';
        X_end_image = (P * X_end.').';
        d = d ./ repmat(estimated_length, 1, 3); % Normalize
        d_image = (P * [d 0].').';
        
        image_line = cross(d_image, X_tip_image);
        u = uFromImageLine(image_line);
        [X_tip_image, X_end_image, tangent_3D] = parametersFromMidline(image_line);
        X_tip_basis_ray = (P_inv * X_tip_image.').';
        X_end_basis_ray = (P_inv * X_end_image.').';
        
        if verbose
            plotProbeReprojection(...
                I, above, below, lengths, widths, P, d, u, X_tip(1:3),...
                sprintf('Linear solution obtained at iteration %d', i)...
            );
        end
    end

updateSolution(p);

% Iteratively refine
l2Norm_past = Inf;
l2Norm = norm(A * p - b);
if verbose
    fprintf('l2Norm for iteration %d = %g\n', i, l2Norm);
    X_tip %#ok<NOPRT>
    d %#ok<NOPRT>
    u %#ok<NOPRT>
end
while (l2Norm_past > l2Norm) &&...
        ((l2Norm_past - l2Norm) / l2Norm > threshold)
    l2Norm_past = l2Norm;
    
    i = i + 1;
    
    A = rhs(X_tip_basis_ray, X_end_basis_ray, tangent_3D);
    b = lhs(u);
    p = A \ b;
    updateSolution(p);
    
    l2Norm = norm(A * p - b);
    
    if verbose
        fprintf('l2Norm for iteration %d = %g\n', i, l2Norm);
        X_tip %#ok<NOPRT>
        d %#ok<NOPRT>
        u %#ok<NOPRT>
    end
    
end

% The probe must be in-front of the camera
X_tip = X_tip(1:3) / X_tip(4);
depth_X_tip = depthFromCamera(P, X_tip);
if depth_X_tip < 0
    X_tip = P_center(1:3) - (X_tip - P_center(1:3));
    camera_zAxis = camera_zAxis.';
    d = d - 2 * dot(camera_zAxis, d) * camera_zAxis;
end

end


