function [ classifier, likelihood ] = mlDiscreteClassifier( class_distributions, increments, varargin )
% MLDISCRETECLASSIFIER  Maximum likelihood classification from discrete distributions
%
% ## Syntax
% classifier = mlDiscreteClassifier(...
%   class_distributions, increments [, 'periodic']...
% )
% classifier = mlDiscreteClassifier(...
%   class_distributions, increments, background_distribution...
% )
% [ classifier, likelihood ] = mlDiscreteClassifier(____)
%
% ## Description
% classifier = mlDiscreteClassifier(...
%   class_distributions, increments [, 'periodic']...
% )
%   Returns a maximum likelihood classifier assuming a uniform, optionally
%   periodic, background distribution.
%
% classifier = mlDiscreteClassifier(...
%   class_distributions, increments, background_distribution...
% )
%   Returns a maximum likelihood classifier using the given background
%   distribution
%
% [ classifier, likelihood ] = mlDiscreteClassifier(____)
%   Additionally returns the likelihood values corresponding to the
%   classification.
%
% ## Input Arguments
%
% class_distributions -- Class density estimators
%   An array, where the last dimension indexes classes.
%   `class_distributions(...,i)` is an array where each element is the
%   conditional probability of a sample value, at the coordinates
%   represented by its position in the array, given the i-th class.
%
%   If there is one class distribution, as can be inferred from the length
%   of `increments`, or the size of `background_distribution`,
%   `class_distributions` is assumed to be entirely occupied by this class
%   distribution, rather than having a final redundant singleton dimension
%   indexing class distributions.
%
% increments -- Sampling increments
%   `increments(i)` is the spacing between samples represented by adjacent
%   indices along the i-th dimension in `class_distributions`. If there are
%   multiple class distributions, then `increments` has length
%   `ndims(class_distributions) - 1`. If there is one class distribution,
%   then `increments` has length `ndims(class_distributions)`.
%
%   The volume of the sample space between samples in the distributions is
%   needed to construct a uniform background distribution. If
%   `background_distribution` is passed, `increments` can be empty.
%
% background_distribution -- Background density estimator
%   An array with the same dimensions as `class_distributions(...,i)`
%   describing the conditional probability of a sample value, at the
%   coordinates represented by its position in the array, given that the
%   sample does not belong to any of the classes.
%
% If `background_distribution` is omitted, a uniform distribution is used
% in its place. If `'periodic'` is passed instead of
% `background_distribution`, the uniform distribution is calculated by
% assuming that the sampling domain is periodic. In other words, the last
% value in each dimension of `class_distributions(...,i)` represents the
% same sampling coordinate as the first value in the dimension.
%
% ## Output Arguments
%
% classifier -- Maximum likelihood classifier
%   An array with the same dimensions as `class_distributions(...,i)` where
%   each value contains the index of the class giving a sample, at the
%   location corresponding to the subscripts of the value, the highest
%   conditional probability. Indices of zero correspond to the background.
%
%   If there is only one class distribution, `classifier` has the same
%   dimensions as `class_distributions`.
%
% likelihood -- Likelihood values
%   An array with the same dimensions as `classifier`, where each value
%   contains the conditional probability of a sample, at the location
%   corresponding to the subscripts of the value, given by the class to
%   which the sample has been assigned by `classifier`. In other words,
%   `likelihood` contains the maximum conditional probability over all
%   classes used to generate `classifier`.
%
% ## Notes
% - This function will produce a Bayes classifier if the input
%   distributions are pre-multiplied by the probabilities of their
%   respective classes.
%
% See also queryDiscretized1DFunction, hueVariableKernelDensityEstimator, hueGaussianDensityEstimator

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created May 26, 2017

nargoutchk(1, 1);
narginchk(2, 3);

periodic = false;
background_distribution = [];

if ~isempty(varargin)
    if ischar(varargin{1})
        if strcmp(varargin{1}, 'periodic')
            periodic = true;
        else
            error('Unrecognized value of second input argument');
        end
    else
        background_distribution = varargin{1};
    end
end

% Create a uniform background distribution, if necessary
if isempty(background_distribution)
    class_distributions_size = size(class_distributions);
    
    if isempty(increments)
        error('`increments` cannot be empty if no background distribution is passed.')
    elseif length(increments) == (ndims(class_distributions) - 1)
        class_distributions_size = class_distributions_size(1:(end - 1));
    elseif length(increments) ~= ndims(class_distributions)
        error('The length of `increments` must be equal to or one less than the number of dimensions in `class_distributions`.')
    end
    
    if periodic
        n_elements = prod(class_distributions_size - 1);
    else
        n_elements = prod(class_distributions_size);
    end
    
    uniform_value = 1 / (n_elements * prod(increments));
    
    if length(class_distributions_size) == 1
        class_distributions_size = [class_distributions_size, 1];
    end

    background_distribution = ones(class_distributions_size) * uniform_value;
end

% Classification
if size(background_distribution, ndims(background_distribution)) == 1
    max_dimension = ndims(background_distribution);
else
    max_dimension = ndims(background_distribution) + 1;
end
all_distributions = cat(max_dimension, background_distribution, class_distributions);
[likelihood, classifier] = max(all_distributions, [], max_dimension);
classifier = classifier - 1;

end

