function [X_tip, d, u] = probeTipAndOrientation( above, below, lengths, widths, P, image_line, threshold, normalize )
%PROBETIPANDORIENTATION Linear estimation of probe tip and orientation
%
% ## Syntax
% [X_tip, d, u] = probeTipAndOrientation(...
%     above, below, lengths, widths, P, image_line, threshold, normalize...
% )
%
% ## Description
% [X_tip, d, u] = probeTipAndOrientation(...
%     above, below, lengths, widths, P, image_line, threshold, normalize...
% )
%   Uses a linear method to calculate the position of the probe tip and the
%   orientation of the probe, given the positions of points on the probe in
%   the image, and given calibration information for the camera and the
%   probe dimensions.
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
%   the image points in `above` and `below`.
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
% ## References
% - R. Hartley and A. Zisserman. Multiple View Geometry in Computer Vision,
%   2nd Edition. Cambridge, UK: Cambridge University Press, 2003.
%
% See also normalizePointsPCA, homography1D

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created March 8, 2017

nargoutchk(3, 3);
narginchk(8, 8);

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
        u_image = u_image(1:2) / repmat(u_image(end), 1, 2);
        y_centerline = -(image_line(1) * u_image(1) + image_line(3)) / image_line(2);
        if y_centerline > u_image(2)
            u = -u;
        end
    end

u = uFromImageLine(image_line);

n = size(above, 1);
above = [above, ones(n, 1)];
below = [below, ones(n, 1)];
allPoints = [above; below];

l = repmat(lengths, 2, 1);

% Find an initial estimate of the probe tip in the image

    function [X_tip_image] = probeTipInImageFromMidline(image_line)
        % Convert the points to 1D coordinates on the estimated line
        image_line_points = closestPointOnLine(repmat(image_line, n, 1), allPoints(:, 1:2));
        tangent = [image_line(2), -image_line(1)];
        tangent = tangent / norm(tangent);
        image_distances = image_line_points - repmat(image_line_points(1, :), n, 1);
        image_distances = dot(image_distances, repmat(tangent, n, 1), 2);
        H = homography1D( l, image_distances, normalize, threshold );
        X_tip_image = (H * [0, 1]).';
    end

[X_tip_image] = probeTipInImageFromMidline(image_line);

% Parameterize the probe tip as a point on the ray through `X_tip_image`,
% with some offset in the direction of `u`.
P_inv = pinv(P); % Pseudoinverse
P_center = null(P); % Camera center
if P_center(end) ~= 0
    P_center = P_center ./ repmat(P_center(end), 1, 4);
end
X_tip_basis_ray = (P_inv * X_tip_image.').';

% From P * (X_tip + l_i * [d; 0] + r_i * [u; 0]) ~ x_i
% So cross(P * (X_tip + l_i * [d; 0] + r_i * [u; 0]), x_i) = 0
% Then substitute X_tip = P_center + lambda1 * X_tip_basis_ray + lambda2 * [u; 0]
%
% Express as a matrix A * [lambda1; lambda2; d(1); d(2); d(3)] = b
x1 = allPoints(:, 1);
x2 = allPoints(:, 2);
r = repmat(widths / 2, 2, 1);
r(n+1:end) = -r(n+1:end); % Account for the opposition between `above` and `below`

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

    function A = rhsFromU(u)
        u1 = u(1);
        u2 = u(2);
        u3 = u(3);        
        A = [
            (x2*(P3_1*P_center1 + P3_2*P_center2 + P3_3*P_center3 + P3_4*P_center4) - P2_1*P_center1 - P2_2*P_center2 - P2_3*P_center3 - P2_4*P_center4),...
            (x2*(P3_1*u1 + P3_2*u2 + P3_3*u3) - P2_2*u2 - P2_3*u3 - P2_1*u1),...
            (P3_1*l*x2 - P2_1*l),...
            (P3_2*l*x2 - P2_2*l),...
            (P3_3*l*x2 - P2_3*l);

            (P1_1*P_center1 - x1*(P3_1*P_center1 + P3_2*P_center2 + P3_3*P_center3 + P3_4*P_center4) + P1_2*P_center2 + P1_3*P_center3 + P1_4*P_center4),...
            (P1_1*u1 + P1_2*u2 + P1_3*u3 - x1*(P3_1*u1 + P3_2*u2 + P3_3*u3)),...
            (P1_1*l - P3_1*l*x1),...
            (P1_2*l - P3_2*l*x1),...
            (P1_3*l - P3_3*l*x1);

            (x1*(P2_1*P_center1 + P2_2*P_center2 + P2_3*P_center3 + P2_4*P_center4) - x2*(P1_1*P_center1 + P1_2*P_center2 + P1_3*P_center3 + P1_4*P_center4)),...
            (x1*(P2_1*u1 + P2_2*u2 + P2_3*u3) - x2*(P1_1*u1 + P1_2*u2 + P1_3*u3)),...
            (P2_1*l*x1 - P1_1*l*x2),...
            (P2_2*l*x1 - P1_2*l*x2),...
            (P2_3*l*x1 - P1_3*l*x2)
            ];
    end

A = rhsFromU(u);

    function b = lhsFromU(u)
        u1 = u(1);
        u2 = u(2);
        u3 = u(3);
        b = [
            x2*(P3_1*(P_center1 + r*u1) + P3_2*(P_center2 + r*u2) + P3_3*(P_center3 + r*u3) + P3_4*P_center4) - P2_1*(P_center1 + r*u1) - P2_2*(P_center2 + r*u2) - P2_3*(P_center3 + r*u3) - P2_4*P_center4;
            P1_1*(P_center1 + r*u1) - x1*(P3_1*(P_center1 + r*u1) + P3_2*(P_center2 + r*u2) + P3_3*(P_center3 + r*u3) + P3_4*P_center4) + P1_2*(P_center2 + r*u2) + P1_3*(P_center3 + r*u3) + P1_4*P_center4;
            x1*(P2_1*(P_center1 + r*u1) + P2_2*(P_center2 + r*u2) + P2_3*(P_center3 + r*u3) + P2_4*P_center4) - x2*(P1_1*(P_center1 + r*u1) + P1_2*(P_center2 + r*u2) + P1_3*(P_center3 + r*u3) + P1_4*P_center4)
        ];
    end

b = lhsFromU(u);

% Minimize L2 norm of (A.p - b)
p = A \ b;

% Update the estimate of `u`
    function updateSolution(p)
        d = p(3:end).';
        d = d ./ repmat(norm(d), 1, 3); % Normalize
        d_image = (P * [d 0].').';
        X_tip = P_center + p(1) * X_tip_basis_ray + p(2) * [u; 0];
        X_tip_image = (P * X_tip.').';
        image_line = cross(d_image, X_tip_image);
        u = uFromImageLine(image_line);
        X_tip_image = probeTipInImageFromMidline(image_line);
        X_tip_basis_ray = (P_inv * X_tip_image.').';
    end

updateSolution(p);

% Iteratively refine
l2Norm_past = Inf;
l2Norm = norm(A * p - b);
while (l2Norm_past > l2Norm) &&...
        ((l2Norm_past - l2Norm) / l2Norm > threshold)
    l2Norm_past = l2Norm;
    
    A = rhsFromU(u);
    b = lhsFromU(u);
    p = A \ b;
    updateSolution(p);
    
    l2Norm = norm(A * p - b);
end

end


