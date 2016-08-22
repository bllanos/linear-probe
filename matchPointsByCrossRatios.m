function [ subject_match_indices ] = matchPointsByCrossRatios( subject, query, subject_gap_cost, query_gap_cost, varargin )
% MATCHPOINTSBYCROSSRATIOS  Align sequences of points on two lines by comparing cross ratios
%
% ## Syntax
% subject_match_indices = matchPointsByCrossRatios(...
%   subject, query, subject_gap_cost, query_gap_cost [, n_samples, verbose]...
% )
%
% ## Description
% subject_match_indices = matchPointsByCrossRatios(...
%   subject, query, subject_gap_cost, query_gap_cost [, n_samples, verbose]...
% )
%   Returns the indices of points in the subject sequence which match
%   points in the query sequence.
%
% ## Input Arguments
%
% subject -- First sequence of collinear of points
%   A vector of length 'm' containing the 1D coordinates of points on a
%   line. The vector should be sorted.
%
% query -- Second sequence of collinear points
%   A vector of length 'n' containing the 1D coordinates of points on the
%   same or a different line from the points in `subject`. The points in
%   `subject` and `query` can be in different coordinate frames. `query`
%   should be sorted.
%
% subject_gap_cost -- Score for gaps in the subject sequence
%   The score to assign to gaps of any length in the subject sequence in
%   the alignment of the subject and query sequences. This value is used
%   indirectly by `swSequenceAlignmentAffine()`.
%
% query_gap_cost -- Score for gaps in the query sequence
%   The score to assign to gaps of any length in the query sequence in
%   the alignment of the subject and query sequences. This value is used
%   indirectly by `swSequenceAlignmentAffine()`.
%
% n_samples - Reduce problem size
%   An integer value from zero to 'n - 1'. If greater than zero, the function
%   will perform fewer computations. The result will be the same provided
%   that the input sequences are composed of points with strongly
%   distinguishable spacings. See below concerning the description of the
%   algorithm for details.
%
%   Defaults to zero if not passed.
%
% verbose -- Debugging flag
%   If true, graphical and console output will be generated for debugging
%   purposes.
%
%   Defaults to false if not passed.
%
% ## Output Arguments
%
% subject_match_indices -- Sequence alignment
%   A vector of length 'n' where `subject_match_indices(i)` is the index of
%   the matched point of `query(i)` in `subject`. If
%   `subject_match_indices(i)` is zero, then `query(i)` is not matched to
%   any point in `subject`.
%
% ## Algorithm
% 
% The two series of points are matched based on their cross ratios, because
% cross ratios (as opposed to positions or length ratios, for instance) are
% invariant to projective transformations:
% 
% First, the sequences of all possible cross ratios (from all possible
% combinations of points) are generated for the two sets of points. These
% sequences of cross ratios are then aligned using
% `swSequenceAlignmentAffine`. A semi-global alignment scheme is used, as
% `query` is assumed to contain only valid points, but possibly be missing
% points, whereas `subject` is assumed to be the sequence of all valid
% points. Consequently, the complete alignment of `query` will be favoured,
% whereas the process will not penalize the extension of `subject` outside
% the aligned region. Lastly, each pair of matched cross ratios in the
% alignment is used to vote for matches between the corresponding points.
% The final point matches are produced using a sequence alignment procedure
% (semi-global alignment with `swSequenceAlignmentAffine`), where the votes
% are matching scores for pairs of points.
%
% If `n_samples` is nonzero, then the sequence of cross ratios generated
% for `query` is constructed as follows, rather than consisting of all
% possible cross ratios: For each point in `query`, a cross ratio is
% computed for the point and a random three other points from `query`. This
% is repeated `n_samples` times, so in total, `n_samples * n` cross ratios are
% used for sequence alignment.
%
% Note that cross ratios are matched with dynamic programming as opposed to
% a closest-pair scheme because a closest-pair scheme would not enforce the
% relative orderings of quadruplets of points on the two lines. The dynamic
% programming algorithm in `swSequenceAlignmentAffine` allows for
% combinations of points to be inserted or deleted, but not transposed.
%
% See also crossRatio, swSequenceAlignmentAffine

% Bernard Llanos
% Spring 2016 research assistantship supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created August 5, 2016

    function [ score ] = f_similarity_cross_ratios_forward(s, q)
        s_val = subject_cross_ratios(s);
        q_val = query_cross_ratios_forward(q);
        score = 1 - abs((q_val - s_val) / s_val);
    end

    function [ score ] = f_similarity_cross_ratios_reverse(s, q)
        s_val = subject_cross_ratios(s);
        q_val = query_cross_ratios_reverse(q);
        score = 1 - abs((q_val - s_val) / s_val);
    end

    function [ score ] = f_similarity_points(s, q)
        score = match_votes(s, q);
    end

nargoutchk(1, 1);
narginchk(4, 6);

n_subject = length(subject);
n_query = length(query);

if ~isempty(varargin)
    n_samples = varargin{1};
    if n_samples < 0 || n_samples > (n_query - 1) || round(n_samples) ~= n_samples
        error('The `n_samples` input argument must be an integer from %d to %d.', 0, n_query - 1)
    end
else
    n_samples = 0;
end
if length(varargin) > 1
    verbose = varargin{2};
else
    verbose = false;
end

% Computing the cross ratio requires 4 points. Disambiguating between
% sequence orientations requires 5 points, because the cross ratio is
% invariant to a reversal of point order.
if n_subject < 5
    error('Insufficient points given on the subject line to compute cross ratios and determine forwards/backwards orientation.')
end
if n_query < 5
    error('Insufficient points given on the query line to compute cross ratios and determine forwards/backwards orientation.')
end

if ~all(subject == sort(subject))
    error('The `subject` input argument should be sorted in ascending order.');
end
if ~all(query == sort(query))
    error('The `query` input argument should be sorted in ascending order.');
end

% Compute all possible cross ratios of the subject and query points
% Note that 'nchoosek' preserves the order of the items being chosen.
subject_combinations = nchoosek(1:n_subject, 4);
n_subject_cross_ratios = size(subject_combinations, 1);
subject_cross_ratios = zeros(n_subject_cross_ratios, 1);
for i = 1:n_subject_cross_ratios
    points = subject(subject_combinations(i, :), 1);
    subject_cross_ratios(i) = crossRatio(points);
end

if n_samples
    
    % Take a random sample of combinations of 4 points
    n_query_cross_ratios = n_samples * n_query;
    triplets = nchoosek(1:(n_query - 1), 3);
    n_triplets = size(triplets, 1);
    sampled_triplets = triplets(randi(n_triplets, n_query_cross_ratios, 1), :);
    query_combinations = zeros(n_query_cross_ratios, 4);
    j = 1;
    for i = 1:n_query
        for k = j:(j + n_samples - 1);
            triplet = sampled_triplets(k, :);
            triplet_filter = triplet < i;
            query_combinations(k, :) = [
                    triplet(triplet_filter),...
                    i,...
                    triplet(~triplet_filter) + 1
                ];
        end
        j = j + n_samples;
    end
    query_combinations = sortrows(query_combinations, 1:4);
else
    
    % Take all combinations of 4 points
    query_combinations = nchoosek(1:n_query, 4);
    n_query_cross_ratios = size(query_combinations, 1);
end

query_cross_ratios_forward = zeros(n_query_cross_ratios, 1);
for i = 1:n_query_cross_ratios
    points = query(query_combinations(i, :), 1);
    query_cross_ratios_forward(i) = crossRatio(points);
end

query_reverse_map = n_query:(-1):1;
query_cross_ratios_reverse = zeros(n_query_cross_ratios, 1);
for i = 1:n_query_cross_ratios
    points = query(query_reverse_map(query_combinations(i, :)), 1);
    query_cross_ratios_reverse(i) = crossRatio(points);
end

if verbose
    % Plot cross ratios
    figure;
    hold on
    plot(subject_cross_ratios, 'g-');
    plot(query_cross_ratios_forward, 'r-');
    hold off
    title('Cross ratios')
    xlabel('Combination index')
    ylabel('Cross ratio')
    legend('Subject sequence cross ratios', 'Query sequence cross ratios')
end

if verbose
    % Plot cross ratios in reverse
    figure;
    hold on
    plot(subject_cross_ratios, 'g-');
    plot(query_cross_ratios_reverse, 'r-');
    hold off
    title('Cross ratios, with query sequence reversed')
    xlabel('Combination index')
    ylabel('Cross ratio')
    legend('Subject sequence cross ratios', 'Reversed query sequence cross ratios')
end

% Match sequences in the given directions with dynamic programming
subject_sequence = 1:n_subject_cross_ratios;
query_sequence = 1:n_query_cross_ratios;
threshold = -Inf;
subject_gap_cost = [subject_gap_cost 0];
query_gap_cost = [query_gap_cost 0];
[ alignment_forward, score_forward ] = swSequenceAlignmentAffine(...
        subject_sequence, query_sequence,...
        @f_similarity_cross_ratios_forward, threshold, subject_gap_cost, query_gap_cost, 'SemiGlobal'...
    );

if verbose
    disp('Forward cross ratio sequence alignment score:');
    disp(score_forward);
end

% Match sequences in the reverse directions with dynamic programming
[ alignment_reverse, score_reverse ] = swSequenceAlignmentAffine(...
        subject_sequence, query_sequence,...
        @f_similarity_cross_ratios_reverse, threshold, subject_gap_cost, query_gap_cost, 'SemiGlobal'...
    );

if verbose
    disp('Reverse cross ratio sequence alignment score:');
    disp(score_reverse);
end

% Vote for point matches using the matched cross ratios
if score_forward > score_reverse
    alignment = alignment_forward;
else
    alignment = alignment_reverse;
end
alignment = alignment(all(alignment, 2), :);

match_votes = zeros(n_subject, n_query);
for i = 1:size(alignment, 1)
    subject_combination = subject_combinations(alignment(i, 1), :);
    query_combination = query_combinations(alignment(i, 2), :);
    for j = 1:4
        match_votes(subject_combination(j), query_combination(j)) = match_votes(subject_combination(j), query_combination(j)) + 1;
    end
end

if verbose
    disp('Subject to query point matching votes:');
    disp('(row = subject point index, column = query point index)');
    disp(match_votes);
end

subject_sequence = 1:n_subject;
query_sequence = 1:n_query;
[ alignment, score ] = swSequenceAlignmentAffine(...
        subject_sequence, query_sequence,...
        @f_similarity_points, threshold, subject_gap_cost, query_gap_cost, 'SemiGlobal'...
    );

if score_forward <= score_reverse
    alignment_filter = (alignment(:, 2) ~= 0);
    alignment(alignment_filter, 2) = query_reverse_map(alignment(alignment_filter, 2));
end

if verbose
    disp('Point sequence alignment:');
    disp(alignment);
    disp('Point sequence alignment score:');
    disp(score);
end

subject_match_indices = zeros(n_query, 1);
for i = 1:n_query
    match_in_subject = alignment(alignment(:, 2) == i, 1);
    if match_in_subject
        subject_match_indices(i) = match_in_subject;
    end
end

end

