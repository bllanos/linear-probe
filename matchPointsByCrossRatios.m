function [ subject_match_indices ] = matchPointsByCrossRatios( subject_line, query_line, subject_gap_cost, query_gap_cost, varargin )
% MATCHPOINTSBYCROSSRATIOS  Align sequences of points on two lines by comparing cross ratios
%
% ## Syntax
% subject_match_indices = matchPointsByCrossRatios(...
%   subject_line, query_line, subject_gap_cost, query_gap_cost [, verbose]...
% )
%
% ## Description
% subject_match_indices = matchPointsByCrossRatios(...
%   subject_line, query_line, subject_gap_cost, query_gap_cost [, verbose]...
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
%   indirectly by `swSequenceAlignment()`.
%
% query_gap_cost -- Score for gaps in the query sequence
%   The score to assign to gaps of any length in the query sequence in
%   the alignment of the subject and query sequences. This value is used
%   indirectly by `swSequenceAlignment()`.
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
% sequences of cross ratios are then aligned using `swSequenceAlignment`. A
% semi-global alignment scheme is used, as `query` is assumed to contain
% only valid points, but possibly be missing points, whereas `subject` is
% assumed to be the sequence of all valid points. Consequently, the
% complete alignment of `query` will be favoured, whereas the process will
% not penalize the extension of `subject` outside the aligned region.
% Lastly, each pair of matched cross ratios in the alignment is used to
% vote for matches between the corresponding points. The final point
% matches are those with the most votes, according to a simple priority
% scheme that prevents many-to-one matching.
%
% Note that cross ratios are matched with dynamic programming as opposed to
% a closest-pair scheme because a closest-pair scheme would not enforce the
% relative orderings of quadruplets of points on the two lines. The dynamic
% programming algorithm in `swSequenceAlignment` allows for combinations of
% points to be inserted or deleted, but not transposed.
%
% See also crossRatio, swSequenceAlignment

% Bernard Llanos
% Spring 2016 research assistantship supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created August 5, 2016

    function [ score ] = f_similarity(s, q)
        s_val = subject_cross_ratios(s);
        q_val = query_cross_ratios(q);
        score = 1 - abs((q_val - s_val) / s_val);
    end

    function [ score ] = f_subject_gap(~)
        score = subject_gap_cost;
    end

    function [ score ] = f_query_gap(~)
        score = query_gap_cost;
    end

nargoutchk(1, 1);
narginchk(4, 5);

if ~isempty(varargin)
    verbose = varargin{1};
else
    verbose = false;
end

n_subject = length(subject_line);
n_query = length(query_line);

% Computing the cross ratio requires 4 points. Disambiguating between
% sequence orientations requires 5 points, because the cross ratio is
% invariant to a reversal of point order.
if n_subject < 5
    error('Insufficient points given on the subject line to compute cross ratios and determine forwards/backwards orientation.')
end
if n_query < 5
    error('Insufficient points given on the query line to compute cross ratios and determine forwards/backwards orientation.')
end

if ~all(subject_line == sort(subject_line))
    error('The `subject` input argument should be sorted in ascending order.');
end
if ~all(query_line == sort(query_line))
    error('The `query` input argument should be sorted in ascending order.');
end

% Compute all possible cross ratios of the subject and query points
% Note that 'nchoosek' preserves the order of the items being chosen.
subject_combinations = nchoosek(1:n_subject, 4);
n_subject_cross_ratios = size(subject_combinations, 1);
subject_cross_ratios = zeros(n_subject_cross_ratios, 1);
for i = 1:n_subject_cross_ratios
    points = subject_line(subject_combinations(i, :), 1);
    subject_cross_ratios(i) = crossRatio(points);
end

query_combinations = nchoosek(1:n_query, 4);
n_query_cross_ratios = size(query_combinations, 1);
query_cross_ratios = zeros(n_query_cross_ratios, 1);
for i = 1:n_query_cross_ratios
    points = query_line(query_combinations(i, :), 1);
    query_cross_ratios(i) = crossRatio(points);
end

if verbose
    % Plot cross ratios
    figure;
    hold on
    plot(subject_cross_ratios, 'g-');
    plot(query_cross_ratios, 'r-');
    hold off
    title('Cross ratios')
    xlabel('Combination index')
    ylabel('Cross ratio')
    legend('Subject sequence cross ratios', 'Query sequence cross ratios')
end

% Match sequences in the given directions with dynamic programming
subject_sequence = 1:n_subject_cross_ratios;
query_sequence_forward = 1:n_query_cross_ratios;
threshold = -Inf;
[ alignment_forward, score_forward ] = swSequenceAlignment(...
        subject_sequence, query_sequence_forward,...
        @f_similarity, threshold, @f_subject_gap, @f_query_gap, 'SemiGlobal'...
    );

% Match sequences in the reverse directions with dynamic programming
query_sequence_reverse = n_query_cross_ratios:-1:1;
threshold = -Inf;
[ alignment_reverse, score_reverse ] = swSequenceAlignment(...
        subject_sequence, query_sequence_reverse,...
        @f_similarity, threshold, @f_subject_gap, @f_query_gap, 'SemiGlobal'...
    );

if verbose
    disp('Forward sequence alignment score:');
    disp(score_forward);
    disp('Reverse sequence alignment score:');
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

subject_match_indices = zeros(n_query, 1);
for i = 1:n_query
    [ max_column_votes, row_indices ] = max(match_votes);
    [ max_vote, col ] = max(max_column_votes);
    row = row_indices(col);
    if max_vote > 0
        subject_match_indices(col) = row;
    end
    match_votes(:, col) = -1;
    match_votes(row, :) = -1;
end

end

