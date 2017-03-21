function [ X_tip, d, u ] = probeTipAndOrientationNonlinear( above, below, lengths, widths, P, X_tip, d, varargin )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

nargoutchk(2, 3);
narginchk(7, 8);

if ~isempty(varargin)
    verbose = varargin{1};
else
    verbose = false;
end

% Initial Guess
d_length = norm(d);
polar = acos(d(3) / d_length);
azimuth = atan2(d(2), d(1));
p0 = [X_tip, polar, azimuth];
if verbose
    disp('Initial guess:')
    disp(p0)
end

n = size(above, 1);
nAll = 2 * n;
allPoints = [above; below];

l = repmat(lengths, 2, 1);
r = repmat(widths / 2, 2, 1); % Take radii, not diameters
r(n+1:end) = -r(n+1:end); % Account for the opposition between `above` and `below`

    function [X_tip, d] = parseSolution( p )
        % p = [ X_tip, polar, azimuth ]
        X_tip = p(1:3);
        polar = p(4);
        azimuth = p(5);
        d = [
            sin(polar) * cos(azimuth),...
            sin(polar) * sin(azimuth),...
            cos(polar)
            ];
    end

    function [ distances ] = costFunction( p )
        [X_tip, d] = parseSolution( p );
        X_tip_image = (P * [X_tip 1].').';
        d_image = (P * [d 0].').';
        image_line = cross(d_image, X_tip_image);
        u = probeNormalFromImageLine(image_line, P);

        reprojected_points = (P * (repmat([X_tip 1], nAll, 1) + l .* repmat([d 0], nAll, 1) + r .* repmat([u 0], nAll, 1)).').';
        reprojected_points = reprojected_points(:, 1:2) ./ repmat(reprojected_points(:, 3), 1, 2);
        
        distances = allPoints - reprojected_points;
        distances = sqrt(dot(distances, distances, 2));
        
        % I could compute the Jacobian to speed up the optimization process
    end

% Invoke optimizer
options = optimoptions('lsqnonlin', 'FunValCheck', 'on');
if verbose
    options.Display = 'iter-detailed';
else
    options.Display = 'none';
end
p = lsqnonlin(@costFunction, p0, [], [], options);

[X_tip, d] = parseSolution( p );
X_tip_image = (P * [X_tip 1].').';
d_image = (P * [d 0].').';
image_line = cross(d_image, X_tip_image);
u = probeNormalFromImageLine(image_line, P);

if verbose
    disp('Final guess:')
    disp(p)
    X_tip
    d
end

end

