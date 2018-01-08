function [ frames ] = readImageGroup(in_directory, wildcard)
% READIMAGEGROUP  Read a sequence of images
%
% ## Syntax
% frames = loadImageGroup(in_directory, wildcard)
%
% ## Description
% frames = loadImageGroup(in_directory, wildcard)
%   Read image files matching the given wildcard, in the given directory,
%   and return them as a 4D array.
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
% frames -- Image sequence
%   A 4D array, where the last dimension is the index of a frame in the
%   sequence of images (in the order listed by `ls()`), and the third
%   dimension indices the colour channels of each image.
%
% ## Notes
% - This function will not work on Windows, as `ls()` has platform-specific
%   behaviour.
%
% ## References
% - Code created in August 2017 for loading RAW images.
% - Code created during Fall 2017 for the CMPUT 551 Mini project on colour
%   classification
%
% See also ls, imread, imwrite

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created January 8, 2018

nargoutchk(1,1)
narginchk(2,2)

% Find all filenames
% `strsplit` will not work if the filenames contain spaces.
names = strsplit(ls(fullfile(in_directory, wildcard)));
n = length(names) - 1; % There is always a terminating newline

% Load the images
I_1 = im2double(imread(names{1}));
sz = size(I_1);
if length(sz) < 3
    sz = [sz, 1];
end
frames = zeros([sz, n]);
frames(:, :, :, 1) = I_1;
for i = 2:n
    I_i = im2double(imread(names{i}));
    frames(:, :, :, i) = I_i;
end

end

