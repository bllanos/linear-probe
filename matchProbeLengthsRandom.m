function [ subject_match_indices ] = matchProbeLengthsRandom(...
    subject_lengths, query_lengths,...
    subject_gap_cost, query_gap_cost,...
    subject_colors, query_colors, direction_threshold, inlier_threshold,...
    varargin...
)
% MATCHPROBELENGTHSRANDOM  Align sequences of points on two lines using colour adjacency relationships and geometric verification
%
% ## Syntax
% subject_match_indices = matchProbeLengths(...
%   subject_lengths, query_lengths, subject_gap_cost, query_gap_cost,...
%   subject_colors, query_colors, direction_threshold, inlier_threshold [, verbose]...
% )
%
% ## Description
% subject_match_indices = matchProbeLengths(...
%   subject_lengths, query_lengths, subject_gap_cost, query_gap_cost,...
%   subject_colors, query_colors, direction_threshold, inlier_threshold [, verbose]...
% )
%   Returns the indices of points in the subject sequence which match
%   points in the query sequence. The matching is based on the colours to
%   the left and right of each point, followed by geometric verification.
%
% ## Input Arguments
%
% subject_lengths -- First sequence of collinear points
%   A column vector of length 'm' containing the 1D coordinates of points
%   on a line. The elements of `subject_lengths` should be sorted.
%
% query_lengths -- Second sequence of collinear points
%   A column vector of length 'n' containing the 1D coordinates of points
%   on the same or a different line from the points in `subject_lengths`.
%   The points in `subject_lengths` and `query_lengths` can be in different
%   coordinate frames. The elements of `query_lengths` should be sorted.
%
% subject_gap_cost -- Score for gaps in the subject sequence
%   The score to assign to gaps of any length in the subject sequence in
%   the alignment of the subject and query sequences. This value is used
%   indirectly by `swSequenceAlignmentAffine()`.
%
%   Scores from zero to one correspond to possible matching scores for
%   individual points in the sequences, whereas negative scores would be
%   worse scores than any matching scores between points.
%
% query_gap_cost -- Score for gaps in the query sequence
%   Analogous to `subject_gap_cost`: The score to assign to gaps of any
%   length in the query sequence.
%
% subject_colors -- Colours for the first sequence of collinear points
%   A matrix of size 'm x 2', where the i-th row describes the colours
%   surrounding the i-th point in `subject_lengths`.
%
%   The first column contains integer labels for the colours preceding the
%   points, whereas the second column contains integer labels for the colours
%   following the points. Use values of zero for unknown colours, which do
%   not match any colour (including themselves), and values of `-1` for
%   arbitrary colours, which half-match any colour.
%
% query_colors -- Colours for the second sequence of collinear points
%   A matrix of size 'n x 2', where the i-th row describes the colours
%   surrounding the i-th point in `query_lengths`. Analogous to
%   `subject_colors`.
%
% direction_threshold -- Threshold for orientation disambiguation
%   As discussed in the 'Algorithm' section below, if the alignment scores
%   for the two possible relative orientations of the subject and query
%   sequences differ by a factor equal to or greater than
%   `direction_threshold`, the lower-scoring orientation is discarded.
%   (Specifically, if the larger score is a factor of `direction_threshold`
%   times or more the smaller score.)
%
% inlier_threshold -- Threshold for selecting inliers
%   As discussed in the 'Algorithm' section below, final matches between
%   the subject and query sequences are filtered to those with colour-based
%   matching scores greater than or equal to `inlier_threshold`.
%
% verbose -- Debugging flag
%   If true, graphical output will be generated for debugging purposes.
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
% Points are first matched using colour information: A table of matching
% scores for all possible pairs of subject and query points is generated.
% The two sequences of points are aligned using
% `swSequenceAlignmentAffine`, with a local alignment scheme.
%
% For a point in the subject sequence, and a point in the query sequence,
% their colour-based matching score is the sum of the scores determined for
% the colours on the left and right. These individual colour scores are
% determined as follows:
% - Same colours: 0.5 points
% - One or both colours is `-1`, and neither colour is zero: 0.25 points
% - One or both colours is zero, or two different colours: 0 points
%
% Consequently, colour-based matching scores are between zero and one.
%
% The alignment is computed in two directions, to account for the two
% possible relative orientations of the point sequences. For each
% direction, multiple optimal alignments (sampled uniformly at random) are
% requested, because there are generally too few colours for unambiguous
% matching, and since there may be errors in colour labelling. Note that
% separate colour-based matching scores are computed for each alignment
% direction.
%
% If one direction produces an alignment score larger than or equal to
% `direction_threshold` times the alignment score of the other direction,
% all of the sample alignments of the other direction are excluded from
% further processing.
%
% In the next stage, geometric verification, triplets of corresponding
% points are selected at random from the optimal alignments. Each triplet
% is used to compute a 1D homography between the subject and query points.
% Matches for the remaining query points are computed by finding the
% closest subject points to the locations of the query points following
% transformation by the homography. The triplet is given a score equal to
% the sum of the colour-based matching scores for all matches.
%
% The triplet selection and scoring process amounts to the RANSAC
% algorithm. When the random sampling of triplets is finished, inlier
% matches are set to those whose colour-based matching scores are greater
% than or equal to `inlier_threshold`. Query points in outlier matches are
% then paired with zeros in `subject_match_indices`.
%
% See also swSequenceAlignmentAffine, matchProbeLengths

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created April 7, 2017

    function [ score ] = f_similarity_forward(s, q)
        score = color_scores(s, q, 1);
    end

    function [ score ] = f_similarity_reverse(s, q)
        score = color_scores(s, q, 2);
    end

    function plotScores(scores, str)
        for index = 1:2
            figure;
            imagesc(scores(:, :, index))
            colorbar
            if index == 1
                title(sprintf('Forward %s', str))
            else
                title(sprintf('Reverse %s', str))
            end
            xlabel('Query sequence index')
            ylabel('Subject sequence index')
        end
    end

nargoutchk(1, 1);
narginchk(8, 9);

% Parse input arguments
if ~isempty(varargin)
    verbose = varargin{1};
else
    verbose = false;
end

n_subject = length(subject_lengths);
subject_sequence = 1:n_subject;
n_query = length(query_lengths);
query_sequence = 1:n_query;

% Computing a 1D homography requires at least 3 points.
if n_subject < 3 || n_query < 3
    verify_matching = false;
    warning('Insufficient points available for geometric verification of alignments.')
else
    verify_matching = true;
end

if ~all(subject_lengths == sort(subject_lengths))
    error('The `subject_lengths` input argument should be sorted in ascending order.');
end
if ~all(query_lengths == sort(query_lengths))
    error('The `query_lengths` input argument should be sorted in ascending order.');
end


subject_colors_grid = repmat(...
    reshape(subject_colors, n_subject, 1, 2),...
    1, n_query, 2 ...
    );
query_colors_grid = cat(3,...
    repmat(reshape(query_colors, 1, n_query, 2), n_subject, 1, 1),...
    repmat(reshape(fliplr(query_colors), 1, n_query, 2), n_subject, 1, 1)...
    );

color_scores = double(subject_colors_grid == query_colors_grid) * 0.5;
color_scores(subject_colors_grid == -1) = 0.25;
color_scores(query_colors_grid == -1) = 0.25;
color_scores(subject_colors_grid == 0) = 0;
color_scores(query_colors_grid == 0) = 0;

% Sum over the scores for colours on the left and right
color_scores = cat(3,...
    sum(color_scores(:, :, 1:2), 3),...
    sum(color_scores(:, :, 3:4), 3)...
    );

if verbose
    plotScores(color_scores, 'colour matching scores');
end

% Match sequences in the given directions with dynamic programming
local_alignment_threshold = 0; % This is kind of a hidden parameter. Perhaps it should be an input argument
subject_gap_cost = [subject_gap_cost 0];
query_gap_cost = [query_gap_cost 0];
n_alignments = n_subject * n_query;
[ alignments_forward, score_forward ] = swSequenceAlignmentAffine(...
        subject_sequence, query_sequence,...
        @f_similarity_forward, local_alignment_threshold,...
        subject_gap_cost, query_gap_cost, 'Local', n_alignments...
    );

if verbose
    disp('Forward sequence alignment score:');
    disp(score_forward);
end

% Match sequences in the reverse directions with dynamic programming
query_sequence_reverse = fliplr(query_sequence);
[ alignments_reverse, score_reverse ] = swSequenceAlignmentAffine(...
        subject_sequence, query_sequence_reverse,...
        @f_similarity_reverse, local_alignment_threshold,...
        subject_gap_cost, query_gap_cost, 'Local', n_alignments...
    );

if verbose
    disp('Reverse sequence alignment score:');
    disp(score_reverse);
end

% Filter out gaps
all_alignments = [alignments_forward, alignments_reverse];
n_all_alignments = 2 * n_alignments;
for i = 1:n_all_alignments
    alignment_filter = all(all_alignments{i}, 2);
    all_alignments{i} = all_alignments{i}(alignment_filter, :);
end

% Correct for direction reversal
for i = (n_alignments + 1):n_all_alignments
    all_alignments{i}(:, 2) = query_sequence_reverse(all_alignments{i}(:, 2));
end

% Visualize the specificity of point matches
if verbose
    votes = zeros(n_subject, n_query, 2);
    for i = 1:n_all_alignments
        alignment = all_alignments{i};
        for j = 1:size(alignment, 1)
            s = alignment(j, 1);
            q = alignment(j, 2);
            if i <= n_alignments
                votes(s, q, 1) = votes(s, q, 1) + 1;
            else
                votes(s, q, 2) = votes(s, q, 2) + 1;
            end
        end
    end
    plotScores(votes, 'matching frequencies');
end

if score_forward >= (direction_threshold * score_reverse)
    all_alignments = all_alignments(1:n_alignments);
    n_all_alignments = n_alignments;
elseif score_reverse >= (direction_threshold * score_forward)
    all_alignments = all_alignments((n_alignments + 1):end);
    n_all_alignments = n_alignments;
end

% TODO Geometric verification
alignment = all_alignments{1};

if verbose
    disp('Sequence alignment:');
    disp(alignment);
end

subject_match_indices = zeros(n_query, 1);
for i = 1:n_query
    match_in_subject = alignment(alignment(:, 2) == i, 1);
    if match_in_subject
        subject_match_indices(i) = match_in_subject;
    end
end

end

