function [ regions_filtered, bw_filtered ] = detectProbeBinaryRegions( regions, bw, radius_adj, varargin )
% DETECTPROBEBINARYREGIONS  Filter binary regions to those which could belong to the probe
%
% ## Syntax
% regions_filtered = detectProbeBinaryRegions(...
%   regions, bw, radius_adj [, verbose]...
% )
% [ regions_filtered, bw_filtered ] = detectProbeBinaryRegions(...
%   regions, bw, radius_adj [, verbose]...
% )
%
% ## Description
% regions_filtered = detectProbeBinaryRegions(...
%   regions, bw, radius_adj [, verbose]...
% )
%   Returns binary regions that satisfy basic probe geometry constraints.
%
% [ regions_filtered, bw_filtered ] = detectProbeBinaryRegions(...
%   regions, bw, radius_adj [, verbose]...
% )
%   Additionally returns binary images corresponding to the binary regions.
%
% ## Input Arguments
%
% regions -- Initial binary regions
%   A structure vector of length n, where the i-th element contains the
%   candidate regions corresponding to the i-th colour. Each element is of
%   the form of the output argument of `bwconncomp`.
%
% bw -- Initial binary images
%   An image_height x image_width x n logical array, where `bw(:,:,i)` is a
%   binary image containing the regions tentatively belonging to the i-th
%   colour class of the probe, contained in `regions(i)`.
%
% radius_adj -- Radius for region proximity test
%   Threshold above which regions are rejected based on their smallest
%   distance to any region belonging to a different colour class.
%
%   This test enforces the requirement that the coloured bands of the probe
%   are adjacent. Note that bands next to occluding objects may be rejected
%   by this test, but this is acceptable, because no usable endpoints can
%   be extracted from such bands for the purposes of probe pose estimation.
%
% verbose -- Debugging flag
%   If true, graphical output will be generated for debugging purposes.
%
%   Defaults to false if not passed.
%
% ## Output Arguments
%
% regions_filtered -- Filtered binary regions
%   A structure vector of length n, where the i-th element contains the
%   regions corresponding to the i-th colour that pass all tests. Each
%   element is of the form of the output argument of `bwconncomp`.
%
% bw_filtered - Filtered binary images
%   An image_height x image_width x n logical array, where
%   `bw_filtered(:,:,i)` is the binary image consisting of the connected
%   components in `regions_filtered(i)`.
%
% See also bwconncomp, extractBinaryRegions

% Bernard Llanos
% Spring 2016 research assistantship supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created August 15, 2016

nargoutchk(1, 2);
narginchk(3, 4);

if ~isempty(varargin)
    verbose = varargin{1};
else
    verbose = false;
end

image_height = size(bw, 1);
image_width = size(bw, 2);
n = size(bw, 3);

regions_filtered = struct('Connectivity', cell(n, 1), 'ImageSize', cell(n, 1),...
    'NumObjects', cell(n, 1), 'PixelIdxList', cell(n, 1));

% Filter regions to those within `radius_adj` of regions for a different
% colour
for i = 1:n
    % Binary image containing all regions for all colours except the i-th
    % colour
    bw_all_but_i = any(bw(:, :, [1:(i-1), (i+1):n]), 3);
    bw_distance = bwdist(bw_all_but_i);
    bw_distance = (bw_distance <= radius_adj);
    if verbose
        figure
        imshow(bw(:, :, i) & bw_distance);
        title(sprintf('Binary image %d AND-ed with thresholded distances to regions for other colours', i))
    end
    bw_filtered_i = false(image_height, image_width);
    regions_i = regions(i);
    region_pixel_lists_i = regions_i.PixelIdxList;
    for j = 1:regions_i.NumObjects
        if any(bw_distance(region_pixel_lists_i{j}))
            bw_filtered_i(region_pixel_lists_i{j}) = true;
        end
    end
    regions_filtered(i) = bwconncomp(bw_filtered_i);
end

if nargout > 1
    bw_filtered = false(image_height, image_width, n);
    for i = 1:n
        regions_filtered_i = regions_filtered(i);
        region_filtered_pixel_lists_i = regions_filtered_i.PixelIdxList;
        bw_filtered_i = false(image_height, image_width);
        for j = 1:regions_filtered_i.NumObjects
            bw_filtered_i(region_filtered_pixel_lists_i{j}) = true;
        end
        bw_filtered(:, :, i) = bw_filtered_i;
        if verbose
            figure
            imshow(bw_filtered_i);
            title(sprintf('Final binary image %d', i))
        end
    end
end

end

