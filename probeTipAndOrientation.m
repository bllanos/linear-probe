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
%   solution.
%
% normalize -- Choice of image point normalization
%   Whether or not to perform normalization of image point coordinates
%   prior to estimating the probe tip and orientation.
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
% See also normalizePointsPCA

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

% Normalize the points
if normalize
    [xNormalized, T] = normalizePointsPCA([above; below]);
    PNormalized = T * P;
else
    xNormalized = [above; below];
    PNormalized = P;
end

% From P * [(X_tip + l_i * d + r_i * u); 1] ~ x_i
% So cross(P * [(X_tip + l_i * d + r_i * u); 1], x_i) = 0
%
% Express as a matrix A * [X_tip(1); X_tip(2); X_tip(3); d(1); d(2); d(3)] = b
x1 = xNormalized(:, 1);
x2 = xNormalized(:, 2);
l = repmat(lengths, 2, 1);
r = repmat(widths / 2, 2, 1);
r(n+1:end) = -r(n+1:end); % Account for the opposition between `above` and `below`

A = [
    (PNormalized(3,1)*x2 - PNormalized(2,1)), (PNormalized(3,2)*x2 - PNormalized(2,2)), (PNormalized(3,3)*x2 - PNormalized(2,3)), (PNormalized(3,1)*l.*x2 - PNormalized(2,1)*l), (PNormalized(3,2)*l.*x2 - PNormalized(2,2)*l), (PNormalized(3,3)*l.*x2 - PNormalized(2,3)*l);
    (PNormalized(1,1) - PNormalized(3,1)*x1), (PNormalized(1,2) - PNormalized(3,2)*x1), (PNormalized(1,3) - PNormalized(3,3)*x1), (PNormalized(1,1)*l - PNormalized(3,1)*l.*x1), (PNormalized(1,2)*l - PNormalized(3,2)*l.*x1), (PNormalized(1,3)*l - PNormalized(3,3)*l.*x1);
    (PNormalized(2,1)*x1 - PNormalized(1,1)*x2), (PNormalized(2,2)*x1 - PNormalized(1,2)*x2), (PNormalized(2,3)*x1 - PNormalized(1,3)*x2), (PNormalized(2,1)*l.*x1 - PNormalized(1,1)*l.*x2), (PNormalized(2,2)*l.*x1 - PNormalized(1,2)*l.*x2), (PNormalized(2,3)*l.*x1 - PNormalized(1,3)*l.*x2)
    ];

    function b = lhsFromU(u)
        b = [
            x2.*(PNormalized(3,4) + PNormalized(3,1)*r*u(1) + PNormalized(3,2)*r*u(2) + PNormalized(3,3)*r*u(3)) - PNormalized(2,4) - PNormalized(2,1)*r*u(1) - PNormalized(2,2)*r*u(2) - PNormalized(2,3)*r*u(3);
            PNormalized(1,4) - x1.*(PNormalized(3,4) + PNormalized(3,1)*r*u(1) + PNormalized(3,2)*r*u(2) + PNormalized(3,3)*r*u(3)) + PNormalized(1,1)*r*u(1) + PNormalized(1,2)*r*u(2) + PNormalized(1,3)*r*u(3);
            x1.*(PNormalized(2,4) + PNormalized(2,1)*r*u(1) + PNormalized(2,2)*r*u(2) + PNormalized(2,3)*r*u(3)) - x2.*(PNormalized(1,4) + PNormalized(1,1)*r*u(1) + PNormalized(1,2)*r*u(2) + PNormalized(1,3)*r*u(3))
        ];
    end

b = lhsFromU(u);

% Minimize L2 norm of (A.p - b)
p = A \ b;

% Update the estimate of `u`
    function updateSolution(p)
        d = p(4:end).';
        d = d ./ repmat(norm(d), 1, 3); % Normalize
        d_image = (P * [d 0].').';
        X_tip = p(1:3).';
        X_tip_image = (P * [X_tip 1].').';
        image_line = cross(d_image, X_tip_image);
        u = uFromImageLine(image_line);
    end

updateSolution(p);

% Iteratively refine
l2Norm_past = Inf;
l2Norm = norm(A * p - b);
while (l2Norm_past > l2Norm) &&...
        ((l2Norm_past - l2Norm) / l2Norm > threshold)
    l2Norm_past = l2Norm;
    
    b = lhsFromU(u);
    p = A \ b;
    updateSolution(p);
    
    l2Norm = norm(A * p - b);
end

end


