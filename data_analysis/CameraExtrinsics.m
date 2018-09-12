%% Computation of camera pose from a new image containing a chequerboard
%
% ## Usage
%   Modify paths to input data in the first code section below, then run.
%
% ## Input
%
% ### Camera calibration
% A '.mat' file containing a 'cameraParams' variable. The variable is a
% structure of class 'cameraParameters' describing the camera, such as
% output by the MATLAB Camera Calibrator app.
%
% ### Image of a chequerboard
% An image of a chequerboard. The camera's pose is to be estimated relative
% to the chequerboard.
%
% ## Output
%
% ### Camera pose
% A '.mat ' file containing the following variables will be saved to a
% location chosen interactively by the user:
% - 'R': The transpose of the `rotationMatrix` output argument of the
%   MATLAB 'extrinsics()' function.
% - 't': The transpose of the `translationVector` output argument of the
%   MATLAB 'extrinsics()' function.
%
% 'R' and 't' represent the Euclidean transformation which converts points
% in world coordinates to points in camera coordinates. The origin in world
% coordinates is the top left corner of the chequerboard.
%
% Additionally, the file contains the values of all parameters in the first
% section of the script below, for reference. (Specifically, those listed
% in `parameters_list`, which should be updated if the set of parameters is
% changed.)
%
% ## References
% - MATLAB Help, "extrinsics", example: 'ComputeExtrinsicsExample.m'
%   (https://www.mathworks.com/help/vision/ref/extrinsics.html)

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created March 14, 2018

%% Input data and parameters

% List of parameters to save with results
parameters_list = {
        'I_filename',...
        'camera_params_filename'...
    };

% Camera calibration
camera_params_filename = '';

% Image of chequerboard
I_filename = '';

% Chequerboard square size
square_size = 106.95; % Must be in the same units as used for camera calibration

%% Load input data

load(camera_params_filename, 'cameraParams');
if ~exist('cameraParams', 'var')
    error('No camera calibration (`cameraParams`) found in ''%s''.', camera_params_filename)
end

I = imread(I_filename);

%% Computations

% Undistort image.
im = undistortImage(I, cameraParams);

% Find reference object in new image.
[image_points, board_size] = detectCheckerboardPoints(im);

% Generate the world coordinates of the checkerboard corners in the
% pattern-centric coordinate system, with the upper-left corner at (0,0).
world_points = generateCheckerboardPoints(board_size, square_size);

% Compute new extrinsics.
[rotationMatrix, translationVector] = extrinsics(image_points, world_points, cameraParams);

R = rotationMatrix.';
t = translationVector.';

%% Save results to a file
save_variables_list = [ parameters_list, {...
        'R',...
        't'...
    } ];
uisave(save_variables_list,'cameraExtrinsics')
