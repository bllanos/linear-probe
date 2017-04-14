function [ axes ] = pcaAxes2D( coeff, mu )
%PCAAXES2D PCA vectors in the original coordinate space of 2D data
%
% ## Syntax
% axes = pcaAxes2D( coeff, mu )
%
% ## Description
% axes = pcaAxes2D( coeff, mu )
%   Returns the homogenous representations of two-dimensional lines
%   corresponding to the PCA basis vectors, in the space of the original
%   data.
%
% ## Input Arguments
%
% coeff -- Principal component coefficients
%   An 2 x 2 array, containing the principal component coefficients for
%   two-dimensional data, as output by MATLAB's `pca` function.
%
% mu -- Estimated means
%   An 1 x 2 vector, containing the estimated means of two-dimensional
%   data, as output by MATLAB's `pca` function.
%
% ## Output Arguments
%
% axes -- Lines representing PCA components
%   An 2 x 3 array, where `axes(i, :)` is the parameters of a line
%   representing the i-th PCA component in the same frame of reference as
%   the original data. `axes(i, :)` is the line described by the equation
%   `axes(i, 1) * x + axes(i, 2) * y + axes(i, 3) = 0`.
%
%   For convenience, lines have been normalized as described here:
%   http://homepages.inf.ed.ac.uk/rbf/CVonline/LOCAL_COPIES/BEARDSLEY/node2.html
%   ("Manipulating Points and Lines," by Bob Fisher,
%    Fri Nov 7 12:08:26 GMT 1997)
%
% See also pca

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created April 9, 2017

dx = coeff(1, :).';
dy = coeff(2, :).';
c = ones(2, 1);
a = c ./ ((dx * mu(2) ./ dy) - mu(1));
b = -a .* dx ./ dy;
% Normalize
scale = sqrt(a .^ 2 + b .^ 2);
axes = [a, b, c] ./ repmat(scale, 1, 3);

end

