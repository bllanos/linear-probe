function [u] = probeNormalFromImageLine(image_line, P)
    % The plane back-projected through the line of the probe
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
