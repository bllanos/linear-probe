function [ frames ] = readVideo( filename, n )
%READVIDEO Read a video file
%
% Input
% filename - The path to the video file
% n - The number of frames to read. (If empty, all frames are read.)
%
% Output
% frames - A 4D array, where the last dimension is the index of a frame
%          in the sequence of frames from the video
%
% CMPUT 615 Lab 1: Motion Estimation and Tracking
% Bernard Llanos, ID: 1505236
% January 24, 2017

v = VideoReader(filename);

k = 1;
frame1 = readFrame(v);
frames = zeros(size(frame1,1), size(frame1,2), size(frame1,3), n);
frames(:,:,:,1) = frame1;
while hasFrame(v) && (isempty(n) || k <= n)
    frames(:,:,:,k) = readFrame(v);
    k = k + 1;
end
frames = frames(:,:,:,1:(k-1));

end

