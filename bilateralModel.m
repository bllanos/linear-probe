function [ model, lengths, axes, model_px, transform ] = bilateralModel( points, point_alignment_outlier_threshold, distinguish_tips )
% BILATERALMODEL  Interpret points as arranged symmetrically around a line segment
%
% ## Syntax
% model = bilateralModel( points, point_alignment_outlier_threshold, distinguish_tips )
% [ model, lengths ] = bilateralModel(...
%   points,...
%   point_alignment_outlier_threshold,...
%   distinguish_tips...
% )
% [ model, lengths, axes ] = bilateralModel(...
%   points,...
%   point_alignment_outlier_threshold,...
%   distinguish_tips...
% )
% [ model, lengths, axes, model_px ] = bilateralModel(...
%   points,...
%   point_alignment_outlier_threshold,...
%   distinguish_tips...
% )
% [ model, lengths, axes, model_px, transform ] = bilateralModel(...
%   points,...
%   point_alignment_outlier_threshold,...
%   distinguish_tips...
% )
%
% ## Description
% model = bilateralModel(...
%   points,...
%   point_alignment_outlier_threshold,...
%   distinguish_tips...
% )
%   Returns the input points in the coordinates of their PCA components,
%   organized into pairs where possible.
%
% [ model, lengths ] = bilateralModel(...
%   points,...
%   point_alignment_outlier_threshold,...
%   distinguish_tips...
% )
%   Additionally returns a summarization of data in `model` with respect to
%   the first PCA component of the input points.
%
% [ model, lengths, axes ] = bilateralModel(...
%   points,...
%   point_alignment_outlier_threshold,...
%   distinguish_tips...
% )
%   Additionally returns the PCA components of the input points in the form of
%   lines in the original coordinate frame of the points.
%
% [ model, lengths, axes, model_px ] = bilateralModel(...
%   points,...
%   point_alignment_outlier_threshold,...
%   distinguish_tips...
% )
%   Additionally returns the original input points organized in the same
%   way as their equivalents in the space of the PCA components.
%
% [ model, lengths, axes, model_px, transform ] = bilateralModel(...
%   points,...
%   point_alignment_outlier_threshold,...
%   distinguish_tips...
% )
%   Additionally returns the transformation matrix that converts the points
%   from the PCA components space to the reference frame of the original
%   data.
%
% ## Input Arguments
%
% points -- 2D point coordinates
%   An n x 2 array, where 'n' is the number of points, and where the
%   columns store point x and y coordinates, respectively. `points` must
%   contain at least three points for 2D principal components analysis to
%   be possible.
%
% point_alignment_outlier_threshold -- Outlier threshold in standard deviations
%   The number of standard deviations from the mean separation of the first
%   coordinates of points paired based on their first coordinates in PCA
%   component space at which a pair is considered to be a false match.
%
% distinguish_tips -- Endpoint identification flag
%   If `distinguish_tips` is true, points at extremal positions along the
%   line corresponding to the first PCA component that have no matches
%   across the line will be returned separately as the tips of the dataset.
%
% ## Output Arguments
%
% model -- Bilaterally-symmetric representation of points
%   A structure with the following fields:
%   - above: Points located above the first axis in the PCA component
%       space, that are paired with points below the first axis based on
%       the similarity of their first coordinates. `above(i)` is paired
%       with `below(i)`.
%   - below: Points located below the first axis in the PCA component
%       space, that are paired with points above the first axis based on
%       the similarity of their first coordinates. `above(i)` is paired
%       with `below(i)`.
%   - unmatched: Points located above or below the first axis in the PCA
%       component space for which no matches were found based on similar
%       first coordinate values. In detail, `unmatched` is the union of
%       the following:
%       - Points whose closest neighbours (by first coordinate in PCA component
%         space) did not select them as their closest neighbours in return.
%       - Points whose closest neighbours (by first coordinate in PCA component
%         space) are at a distance of more than
%         `point_alignment_outlier_threshold` times the standard deviation
%         of closest neighbour separations from the mean closest neighbour
%         separation.
%       `unmatched` may be an empty array.
%
%   All fields store k x 2 arrays, where k varies, but is equal to the
%   number of points meeting the criteria corresponding to the field, and
%   the columns store the x and y coordinates of points in the coordinate
%   space defined by the PCA components of `points`. `above` has the same
%   dimensions as `below`, because these points are paired. All arrays are
%   sorted by their first column.
%
%   The `distinguish_tips` flag allows for special treatment of extremal
%   points, and may result in additional fields:
%   - If `distinguish_tips` is true, and the first point in `unmatched` has the
%     lowest first coordinate value of all points, then that point is
%     removed from `unmatched` and placed in a field called `head`.
%   - Next, if `distinguish_tips` is true, and the last point in `unmatched` has
%     the highest first coordinate value of all points, then that point is
%     removed from `unmatched` and placed in a field called `tail`.
%
% lengths -- Unification of matched points along the first PCA component line
%   A column vector storing the averaged first PCA component coordinates of the
%   paired points in `model.above` and `model.below`, as well as extremal
%   points, expressed relative to the leftmost point(s).
%
%   Specifically, let 's' equal `model.head(1)`, if `model.head` exists, or
%   otherwise, the average of `model.above(1, 1)` and `model.below(1, 1)`.
%
%   The first element of `lengths` has a value of zero and corresponds to
%   `model.head`, if `model.head` exists, or otherwise, to the pair of
%   `model.above(0, :)` and `model.below(0, :)`. The last element of
%   `lengths` has a value of `model.tail(1) - s` and corresponds to
%   `model.tail`, if `model.tail` exists. Otherwise, the last element of
%   `lengths` corresponds to the pair `model.above(end, :)` and
%   `model.below(end, :)` and is equal to
%   `mean([model.above(end, 1), model.below(end, 1)]) - s`.
%
%   The k-th intermediate element of `lengths` corresponds to the pair
%   `model.above(k, :)` and `model.below(k, :)` and is equal to
%   `mean([model.above(k, 1), model.below(k, 1)]) - s`.
%
% axes -- Lines representing PCA components
%   An 2 x 3 array, where `axes(i, :)` is the parameters of a line
%   representing the i-th PCA component of `points` in the same frame of
%   reference as `points`. `axes(i, :)` is the line described by the
%   equation `axes(i, 1) * x + axes(i, 2) * y + axes(i, 3) = 0`.
%
%   For convenience, lines have been normalized as described here:
%   http://homepages.inf.ed.ac.uk/rbf/CVonline/LOCAL_COPIES/BEARDSLEY/node2.html
%   ("Manipulating Points and Lines," by Bob Fisher,
%    Fri Nov 7 12:08:26 GMT 1997)
%
% model_px -- Bilaterally-symmetric representation of points in pixel space
%   The equivalent of `model`, but the points are in pixel space as opposed
%   to the reference frame defined by their PCA components.
%
%   `model_px` contains the original data points, rather than points
%   obtained using the transformation between the PCA and pixel reference
%   frames. As such, the coordinates of the points are unaffected by
%   numerical error.
%
% transform -- Map back from PCA component space
%   A 3 x 2 array, of the form `[coeff mu; 0 0 1].'`, where `coeff` is the
%   principal components expressed in the reference frame of the original
%   data, and `mu` is the vector of variable means. `transform` can be used
%   to map points in `model` back into their original frame of reference
%   using the expression `[model_points, 1] * transform` (where '1' is a
%   column vector of ones).
%
% ## Notes
% - "Bilateral symmetry" refers to the fact that the function attempts to find
%   pairs of points that are approximately at the same position along the
%   first PCA component line, but which are on opposite sides of the first PCA
%   component line. There is no symmetry in the sense that pairs of points
%   may be at unequal distances from the first PCA component line.
% - Suppose the first PCA component line happens to be centered between pairs of
%   points on either side, and the points represent edge points on an image
%   of a linear object. This does not mean that the PCA component line
%   corresponds to the midline of the object, due to projective distortion.
%
% See also pca

% Bernard Llanos
% Spring 2016 research assistantship supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created August 2, 2016

nargoutchk(1, 5);
narginchk(3, 3);

% Express the points in the coordinate space of their principal components
% using Principal Components Analysis
if nargout > 2
    [ coeff, score, ~, ~, ~, mu ] = pca(points);
    
    % Express the PCA component vectors as lines in the space of the original data
    axes = pcaAxes2D( coeff, mu );
    
    transform = [coeff, mu.'; 0 0 1].';
else
    [~,score] = pca(points);
end

% Keep the points with their PCA component space equivalents, if the
% original points are to be returned as `model_px`
if nargout > 3
    score = [score, points];
end

% Partition points into those above and below the first component
above_line_filter = score(:, 2) > 0;
points_above = score(above_line_filter, :);
points_below = score(~above_line_filter, :);

% Sort points by first coordinate
points_above = sortrows(points_above);
points_below = sortrows(points_below);

% Organize the points into a single row by matching points above and below
% the first component. Unmatched points will be kept separate.
n_points_above = size(points_above, 1);
min_separations_above = zeros(n_points_above, 1);
min_x_separation_indices_above = zeros(n_points_above, 1);
for i = 1:n_points_above
    [min_separations_above(i), min_x_separation_indices_above(i)] =...
        min(abs(points_below(:, 1) - points_above(i, 1)));
end

n_points_below = size(points_below, 1);
min_separations_below = zeros(n_points_below, 1);
min_x_separation_indices_below = zeros(n_points_below, 1);
for i = 1:n_points_below
    [min_separations_below(i), min_x_separation_indices_below(i)] =...
        min(abs(points_above(:, 1) - points_below(i, 1)));
end

% Identify mismatched points based on mutual selection
inliers_above_filter = (min_x_separation_indices_below(min_x_separation_indices_above) == (1:n_points_above).');
inliers_below_filter = (min_x_separation_indices_above(min_x_separation_indices_below) == (1:n_points_below).');

% Identify mismatched points with a statistical analysis
% Reference: MATLAB documentation page on "Inconsistent Data"
% Outlier separation in x-coordinates
inlier_separations = [min_separations_above(inliers_above_filter); min_separations_below(inliers_below_filter)];
sigma_separations = std(inlier_separations);
if sigma_separations > 0
    mu_separations = mean(inlier_separations);

    mu_separations_rep = repmat(mu_separations,n_points_above,1);
    sigma_separations_rep = repmat(sigma_separations,n_points_above,1);
    inliers_above_filter = inliers_above_filter & (abs(min_separations_above - mu_separations_rep) < point_alignment_outlier_threshold * sigma_separations_rep);

    mu_separations_rep = repmat(mu_separations,n_points_below,1);
    sigma_separations_rep = repmat(sigma_separations,n_points_below,1);
    inliers_below_filter = inliers_below_filter & (abs(min_separations_below - mu_separations_rep) < point_alignment_outlier_threshold * sigma_separations_rep);
end

% Organize mismatched points according to their positions along the first
% component axis
unmatched = [
        points_above(~inliers_above_filter, :);
        points_below(~inliers_below_filter, :)
    ];
unmatched = sortrows(unmatched);

% Output the points
model = struct(...
        'above', {points_above(inliers_above_filter, :)},...
        'below', {points_below(inliers_below_filter, :)}...
    );
if distinguish_tips && ~isempty(unmatched)
    if unmatched(1, 1) < model.above(1, 1) && unmatched(1, 1) < model.below(1, 1)
        model.head = unmatched(1, :);
        unmatched = unmatched(2:end, :);
    end
    if ~isempty(unmatched) && unmatched(end, 1) > model.above(end, 1) && unmatched(end, 1) > model.below(end, 1)
        model.tail = unmatched(end, :);
        unmatched = unmatched(1:(end - 1), :);
    end
end
model.unmatched = unmatched;

if nargout > 1
    n_intermediate_lengths = size(model.above, 1);
    start_offset = isfield(model, 'head');
    end_offset = isfield(model, 'tail');
    lengths = zeros(n_intermediate_lengths + start_offset + end_offset, 1);
    if start_offset
        s = model.head(1);
    else
        s = mean([model.above(1, 1), model.below(1, 1)]);
    end
    for i = (start_offset + 1):(n_intermediate_lengths + start_offset)
        lengths(i) = mean([
                model.above(i - start_offset, 1),...
                model.below(i - start_offset, 1)
            ]) - s;
    end
    if end_offset
        lengths(end) = model.tail(1) - s;
    end
end

if nargout > 3
    % Separate points in PCA components space from points in pixel space
    model_px = struct;
    for i = fieldnames(model).'
        field = i{1};
        model_px.(field) = model.(field)(:, 3:4);
        model.(field) = model.(field)(:, 1:2);
    end
end

end