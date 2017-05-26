function [ classifier ] = mlDiscreteClassifier( class_distributions, varargin )
% MLDISCRETECLASSIFIER  Maximum likelihood classification from discrete distributions
%
% ## Syntax
% classifier = mlDiscreteClassifier( class_distributions [, 'single'] )
% classifier = mlDiscreteClassifier( class_distributions, 'periodic' [, 'single'] )
% classifier = mlDiscreteClassifier( class_distributions, background_distribution )
%
% ## Description
% classifier = mlDiscreteClassifier( class_distributions [, 'single'] )
%   Returns a maximum likelihood classifier assuming a uniform background
%   distribution.
%
% classifier = mlDiscreteClassifier( class_distributions, 'periodic' [, 'single'] )
%   Returns a maximum likelihood classifier assuming a uniform, periodic
%   background distribution.
%
% classifier = mlDiscreteClassifier( class_distributions, background_distribution )
%   Returns a maximum likelihood classifier using the given background
%   distribution
%
% ## Input Arguments
%
% class_distributions -- Class density estimators
%   An array, where the last dimension indexes classes.
%   `class_distributions(...,i)` is an array where each element is the
%   conditional probability of a sample value, at the coordinates
%   represented by its position in the array, given the i-th class.
%
%   If the last argument of the function is 'single',
%   `class_distributions` sill be interpreted as a single class
%   distribution. (It is difficult to generate a multidimensional array
%   where the last dimension has size 1. For example, `zeros(10,10,1)`
%   produces a two-dimensional array.)
%
%   Note that 'single' is not passed if `background_distribution` is
%   passed, because it is possible to infer whether `class_distributions`
%   represents a single class distribution, by comparing the dimensions of
%   the two input arguments.
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
%   each value contains the index of the class giving a sample value at the
%   location corresponding to the subscripts of the value the highest
%   conditional probability. Indices of zero correspond to the background.
%
%   If the last argument of the function is 'single', `classifier` has the
%   same dimensions as `class_distributions`.
%
% ## Notes
% - This function will produce a Bayes classifier if the input
%   distributions are multiplied by the probabilities of their respective
%   classes.
%
% See also queryDiscretized1DFunction

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created May 26, 2017

nargoutchk(1, 1);
narginchk(1, 3);

background_distribution = [];
periodic = false;
single = false;

if ~isempty(varargin)
    if ischar(varargin{1})
        if strcmp(varargin{1}, 'periodic')
            periodic = true;
            if length(varargin) > 1
                if strcmp(varargin{1}, 'single')
                    single = true;
                else
                    error('Unrecognized value of third input argument');
                end
            end
        elseif strcmp(varargin{1}, 'single')
            single = true;
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
    if ~single
        class_distributions_size = class_distributions_size(1:(end - 1));
    end
    
    if periodic
        n_elements = prod(class_distributions_size - 1);
    else
        n_elements = prod(class_distributions_size);
    end
    
    if length(class_distributions_size) == 1
        class_distributions_size = [class_distributions_size, 1];
    end

    background_distribution = ones(class_distributions_size) / n_elements;
end

% Classification
if size(background_distribution, ndims(background_distribution)) == 1
    max_dimension = ndims(background_distribution);
else
    max_dimension = ndims(background_distribution) + 1;
end
all_distributions = cat(max_dimension, background_distribution, class_distributions);
[~, classifier] = max(all_distributions, [], max_dimension);
classifier = classifier - 1;

end

