function [H] = homography1D( x1, x2, normalize, varargin )
%HOMOGRAPHY1D Linear estimation of a 1D projective homography
%
% ## Syntax
% [H] = homography1D( x1, x2, normalize [, threshold] )
%
% ## Description
% [H] = homography1D( x1, x2, normalize [, threshold] )
%   Computes a 1D projective homography relating `x1` and `x2` by
%   minimizing algebraic error.
%
% ## Input Arguments
%
% x1 -- First set of points
%   A vector of length n containing the 1D coordinates of points on a line.
%
% x2 -- Second set of points
%   A vector of length n containing the 1D coordinates of points on a
%   second line.
%
% normalize -- Choice of image point normalization
%   Whether or not to perform normalization of point coordinates prior to
%   estimating the homography.
%
% threshold -- Convergence threshold for L1 norm minimization
%   If passed, minimize the L1 norm of the algebraic error, as opposed to
%   the L2 norm of the algebraic error. Minimization is an iterative
%   procedure: When the relative change in the residual error between
%   iterations of the linear solution drops below this threshold, return
%   the current solution.
%
% ## Output Arguments
%
% H -- 1D projective homography
%   A homography such that `x2 ~ (H * x1.').'`, where `x1` and `x2` have
%   been expanded with an additional columns of ones. `H` is a 2 x 2 array.
%
% ## References
% - R. Hartley and A. Zisserman. Multiple View Geometry in Computer Vision,
%   2nd Edition. Cambridge, UK: Cambridge University Press, 2003.

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created March 9, 2017

nargoutchk(1,1);
narginchk(4,4);

n = length(x1);

% Convert to homogenous coordinates
if size(x1, 1) == 1
    x1 = [x1.' ones(n, 1)];
else
    x1 = [x1 ones(n, 1)];
end
if size(x2, 1) == 1
    x2 = [x2.' ones(n, 1)];
else
    x2 = [x2 ones(n, 1)];
end

% Normalize the points
if normalize
    [x1Normalized, T1] = normalizePointsPCA(x1);
    [x2Normalized, T2] = normalizePointsPCA(x2);
else
    x1Normalized = x1;
    x2Normalized = x2;
end

% From lambda_i * x2_i = H * x1_i
A = [
    x1Normalized(:, 1), ones(n, 1), 0, 0, diag(-x2Normalized(:, 1));
    0, 0, x1Normalized(:, 1), ones(n, 1), -eye(n)
    ];

% Minimize L2 norm of A.h
% This provides the initial estimate for L1 minimization as well.
[~,~,V] = svd(A);
h = V(:, end);
    
if ~isempty(varargin)
    % Minimize L1 norm of A.h
    % See 'l1DLT.pdf' provided with the lab instructions
    eta = abs(A * h);
    eta(eta == 0) = 1; % Avoid zero weights
    l1Norm_past = Inf;
    l1Norm = sum(eta);
    threshold = varargin{1};
    while (l1Norm_past > l1Norm) &&...
            ((l1Norm_past - l1Norm) / l1Norm > threshold)
        l1Norm_past = l1Norm;
        etaA = diag((sqrt(eta)).^(-1)) * A;
        [~,~,V] = svd(etaA);
        h = V(:, end);
        eta_new = abs(A * h);
        % Avoid zero weights
        eta(eta_new ~= 0) = eta_new(eta_new ~= 0);
        l1Norm = sum(eta);
        % disp(l1Norm);
    end
end

H = reshape(h(1:4), 2, 2)';

% Undo normalization
if normalize
    H = T2 \ H * T1;
end
end