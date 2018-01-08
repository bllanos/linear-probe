function [px_mean, px_std, sz] = loadImageGroup(in_directory, wildcard)
% LOADIMAGEGROUP  Read and average a group of images
%
% ## Syntax
% [px_mean, px_std, sz] = loadImageGroup(in_directory, wildcard)
%
% ## Description
% [px_mean, px_std, sz] = loadImageGroup(in_directory, wildcard)
%   Read image files matching the given wildcard, in the given directory,
%   and return their average and standard deviation, and pixel location
%   information.
%
% ## Input Arguments
%
% in_directory -- Input directory
%   A character vector containing the path of the directory containing
%   image files to load.
%
% wildcard -- Input filename wildcard
%   A wildcard expression for `ls()`. `wildcard` determines which files in
%   the input directory will be read. All image files must have the same
%   pixel dimensions.
%
%   Filenames must not contain spaces.
%
% ## Output Arguments
%
% px_mean -- Mean pixel values
%   An n x 3 array, where 'n' is the number of pixels in an image,
%   containing the averaged pixels in the images. Each row is one pixel,
%   with the columns containing the Red, Green, and Blue colour channels.
%
% px_std -- Pixel standard deviations
%   An n x 3 array, where 'n' is the number of pixels in an image,
%   containing the standard deviations of pixel values in the images. Each
%   row corresponds to one pixel location, with the columns containing the
%   standard deviations of the Red, Green, and Blue colour channels, for
%   pixels at that location.
%
% sz -- Image dimensions
%   A two-element vector containing the height and width of the images.
%   `sz` can be used to reshape the `px_mean` output argument into an image
%   for display.
%
% ## Notes
%
% ## References
% - Code created in August 2017 for loading RAW images.
%
% See also ls, imread, imwrite

% Created for: CMPUT 551 Mini project: Colour classification
% Fall 2017
% Bernard Llanos, ID 1505236
% Department of Computing Science, University of Alberta

nargoutchk(0,3)
narginchk(2,2)

% Find all filenames
% Caution: This will not work on Windows, as `ls()` has platform-specific behaviour.
% Also, `strsplit` will not work if the filenames contain spaces.
names = strsplit(ls(fullfile(in_directory, wildcard)));
n = length(names) - 1; % There is always a terminating newline

% Load the images
I_1 = im2double(imread(names{1}));
sz = [size(I_1, 1), size(I_1, 2)];
n_px = prod(sz);
I = zeros(n_px, 3, n);
I(:, :, 1) = reshape(I_1, n_px, 3);
for i = 2:n
    I_i = im2double(imread(names{i}));
    I(:, :, i) = reshape(I_i, n_px, 3);
end

px_mean = mean(I, 3);
px_std = std(I, 0, 3);

end

