function [ points, search_windows, corners ] = extractInterestPoints( I_annotation, varargin )
% EXTRACTINTERESTPOINTS  Read and refine user-marked corner points
%
% ## Syntax
% points = extractInterestPoints( I_annotation )
% points = extractInterestPoints( I_annotation, I, search_width )
% [ points, search_windows ] = extractInterestPoints( I_annotation, I, search_width )
% [ points, search_windows, corners ] = extractInterestPoints( I_annotation, I, search_width )
%
% ## Description
% points = extractInterestPoints( I_annotation )
%   Returns the centers of nonzero areas in the input image.
%
% points = extractInterestPoints( I_annotation, I, search_width )
%   Returns corner feature locations within search windows of the centers
%   of nonzero areas in the input annotation image.
%
% [ points, search_windows ] = extractInterestPoints( I_annotation, I, search_width )
%   Additionally returns the search windows used to find corner features.
%
% [ points, search_windows, corners ] = extractInterestPoints( I_annotation, I, search_width )
%   Additionally returns all corner features within the search windows.
%
% ## Input Arguments
%
% I_annotation -- Annotation image
%   An h x w array, containing mostly zero entries and a few nonzero blobs.
%   Conceptually, `I_annotation` stores user-marked interest points.
%
% I -- Base image
%   An h x w x (1|3) array representing the image corresponding to the
%   annotations in `I_annotation`.
%
% search_width -- Search window half width
%   An positive integer specifying how far from each center of an
%   annotation to search for corner feature points. A square search window
%   of at most size `2 x search_width + 1` will be constructed around each
%   center pixel. Search windows will be constrained in size such that the
%   search windows of multiple annotation centers do not overlap.
%
% ## Output Arguments
%
% points -- Adjusted annotation points
%   An n x 2 array, where 'n' is the number of distinct nonzero regions in
%   `I_annotation`, and where the columns store pixel x and y coordinates,
%   respectively.
%
%   The function computes an initial value for `points` by applying MATLAB's
%   morphological 'shrink' operation to `I_annotation` and extracting the
%   coordinates of the resulting nonzero pixels. If `I` and `search_width`
%   are passed, then the locations in `points` are adjusted to be the
%   strongest corner features within the surrounding search windows. If a
%   search window around a location in `points` does not contain any
%   significant features, the location is returned as is.
%
% search_windows -- Corner feature search regions
%   An n x 4 array, where `search_windows(i, :)` is the rectangular region
%   of interest used to locate `points(i, :)`. `search_windows(i, :)` is of
%   the form `[x y width height]`.
%
% corners -- Detected corner features
%   An n x 1 cell vector, where `corners{i}` contains all corner features
%   detected in `search_windows(i, :)`, the strongest of which is
%   `points(i, :)`. `corners{i}` is a `cornerPoints` object.
%
% ## Notes
% - Corner features are currently detected with MATLAB's implementation of
%   minimum eigenvalue algorithm, `detectMinEigenFeatures`, on a greyscale
%   version of `I` obtained with `rgb2gray`.
%
% See also bwmorph, detectMinEigenFeatures, cornerPoints, rgb2gray

% Bernard Llanos
% Spring 2016 research assistantship supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created August 2, 2016

nargoutchk(1, 3);
narginchk(1, 3);
if nargin == 2
    error('`I` and `search_width` must either both be passed or both omitted.');
end
if nargout > 1 && nargin ~= 3
    error('`Search windows cannot be output unless `I` and `search_width` are passed.');
end

I_annotation = logical(I_annotation);
[annotations_y, annotations_x] = find(bwmorph(I_annotation, 'shrink', Inf));
points = [annotations_x, annotations_y];

if nargin > 1
    I = varargin{1};
    if size(I, 3) > 1
        I_grey = rgb2gray(I);
    else
        I_grey = I;
    end
    image_width = size(I, 2);
    image_height = size(I, 1);
    
    search_width = varargin{2};
    if search_width < 0 || round(search_width) ~= search_width
        error('`search_width` must be an integer greater than zero.');
    end
    n_annotations = length(annotations_x);
    
    corners = cell(n_annotations, 1);
    search_windows = zeros(n_annotations, 4);
    for i = 1:n_annotations
        x_separation = annotations_x([1:(i-1), (i+1):end]) - annotations_x(i);
        y_separation = annotations_y([1:(i-1), (i+1):end]) - annotations_y(i);
        distances = sqrt(x_separation .^ 2 + y_separation .^ 2);
        search_width_i = min(floor((min(distances) - 1) / 2), search_width);
        search_width_i = max(0, search_width_i);
        max_width = search_width_i * 2 + 1;
        x_min = max(1, annotations_x(i) - search_width_i);
        width_x = min(max_width, max(1, image_width - x_min));
        y_min = max(1, annotations_y(i) - search_width_i);
        width_y = min(max_width, max(1, image_height - y_min));
        
        search_windows(i, :) = [x_min, y_min, width_x, width_y];
        corners{i} = detectMinEigenFeatures(I_grey, 'ROI', search_windows(i, :));
        
        if ~isempty(corners{i})
            points(i, :) = corners{i}.selectStrongest(1).Location;
        end
    end
end



end