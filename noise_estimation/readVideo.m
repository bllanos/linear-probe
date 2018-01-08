function [ frames ] = readVideo( filename, varargin )
% READVIDEO  Read a sequence of video frames
%
% ## Syntax
% frames = readVideo(filename)
% frames = readVideo(filename, start)
% frames = readVideo(filename, start, n)
%
% ## Description
% frames = readVideo(filename)
%   Reads frames from the given video file and returns them as a 4D array.
% frames = readVideo(filename, start)
%   Reads frames from the given video file, starting from the given offset,
%   and returns them as a 4D array.
% frames = readVideo(filename, start, n)
%   Reads the given number of frames from the given video file, starting
%   from the given offset, and returns them as a 4D array.
%
% ## Input Arguments
%
% filename -- Input video
%   A character vector containing the path of the video file.
%
% start -- First frame
%   The time in seconds of the first frame to read. Defaults to the start
%   of the video if zero, or if not passed.
%
% n -- Number of frames
%   The maximum number of frames to read. Defaults to the number of frames
%   in the rest of the video after `start` if zero, or if not passed.
%
% ## Output Arguments
%
% frames -- Image sequence
%   A 4D array, where the last dimension is the index of a frame in the
%   sequence of frames from the video, and the third dimension indices the
%   colour channels of each image.
%
% ## Notes
%
% ## References
% - Code created during Winter 2017 for CMPUT 615 Lab 1 on Motion
%   Estimation and Tracking
%
% See also VideoReader

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created January 8, 2018

nargoutchk(1,1)
narginchk(1,3)

v = VideoReader(filename);

stop = v.Duration;
frame_rate = v.FrameRate;
if frame_rate <= 0
    v = VideoReader(filename); % The object must be recreated when querying 'NumberOfFrames'
    frame_rate = v.NumberOfFrames / stop; %#ok<VIDREAD>
    v = VideoReader(filename); % The object must be recreated when querying 'NumberOfFrames'
end

if nargin > 1
    if varargin{1} > 0 && varargin{1} < v.Duration
        v.CurrentTime = varargin{1};
    elseif varargin{1} >= v.Duration
        error('Cannot start reading past the end of the video.')
    elseif varargin{1} < 0
        error('Cannot start reading before the start of the video.')
    end
end

if nargin > 2
    if varargin{2} > 0
        n = varargin{2};
        stop = v.CurrentTime + (n / frame_rate);
        if stop > v.Duration
            stop = v.Duration;
        end
    elseif varargin{2} < 0
        error('Cannot read a negative number of frames.')
    end
end
n = ceil((stop - v.CurrentTime) * frame_rate);

frame1 = im2double(readFrame(v));
frames = zeros(size(frame1,1), size(frame1,2), size(frame1,3), n);
frames(:,:,:,1) = frame1;
k = 2;
while hasFrame(v) && k <= n
    frames(:,:,:,k) = im2double(readFrame(v));
    k = k + 1;
end
if k <= n
    frames = frames(:, :, :, 1:(k - 1));
end

end

