function [ val ] = queryDiscretized1DFunction( x, f, inc, varargin )
% QUERYDISCRETIZED1DFUNCTION  Evaluate a sampled function of one argument
%
% ## Syntax
% val = queryDiscretized1DFunction( x, f, inc [, offset ] )
%
% ## Description
% val = queryDiscretized1DFunction( x, f, inc [, offset ] )
%   Returns the values of the discretized function that most closely
%   correspond to the query values.
%
% ## Input Arguments
%
% x -- Query values
%   An array (of any dimensions) containing the values at which to sample the
%   function discretization.
%
% f -- Evaluated function
%   A vector of length m representing the evaluation of a function at
%   equally-spaced 1D points.
%
% inc -- Sample point spacing
%   A scalar representing the (signed) constant spacing between the points
%   at which the function was evaluated. For instance, the k-th element of
%   `f` is equal to the evaluation of the function at a point `inc` units
%   in the positive direction from the point at which the function was
%   evaluated to produce the (k-1)-th element of `f`.
%
% offset -- Sampling axis origin
%   The value at which the function was sampled to produce the first
%   element of `f`.
%
%   Defaults to zero if not passed.
%
% ## Output Arguments
%
% val -- Approximate function evaluation
%   An array with the same dimensions as `x` containing the values from `f`
%   corresponding to the closest points to the points in `x` at which the
%   function was evaluated.
%
% ## Notes
% - For simplicity, this is a nearest-neighbour approach (i.e. zero-order
%   interpolation of the discretized function).

% Bernard Llanos
% Spring 2016 research assistantship supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created August 9, 2016

nargoutchk(1, 1);
narginchk(3, 4);

if ~isempty(varargin)
    offset = varargin{1};
else
    offset = 0;
end

x_vector = x(:);
indices = floor(((x_vector - offset) / inc) + 0.5) + 1;
if any(indices < 0 | indices > length(f))
    error('Some of the query values in input `x` are out of the domain sampled to construct input `f`.')
end
val = f(indices);
val = reshape(val, size(x));

end

