function [ subject_match_indices ] = matchProbeLengths( subject_lengths, query_lengths, subject_gap_cost, query_gap_cost, varargin )
% MATCHPROBELENGTHS  Align sequences of points on two lines using cross ratios and colour adjacency relationships
%
% ## Syntax
% subject_match_indices = matchProbeLengths(...
%   subject_lengths, query_lengths, subject_gap_cost, query_gap_cost...
%   [, verbose]...
% )
%
% subject_match_indices = matchProbeLengths(...
%   subject_lengths, query_lengths, subject_gap_cost, query_gap_cost,...
%   subject_colors, query_colors, color_weight...
%   [, verbose]...
% )
%
% ## Description
% subject_match_indices = matchProbeLengths(...
%   subject_lengths, query_lengths, subject_gap_cost, query_gap_cost...
%   [, verbose]...
% )
%   Returns the indices of points in the subject sequence which match
%   points in the query sequence. The matching is based on point cross
%   ratios.
%
% subject_match_indices = matchProbeLengths(...
%   subject_lengths, query_lengths, subject_gap_cost, query_gap_cost,...
%   subject_colors, query_colors, color_weight...
%   [, verbose]...
% )
%   Returns the indices of points in the subject sequence which match
%   points in the query sequence. The matching is based on point cross
%   ratios and on the colours to the left and right of each point.
%
% ## Input Arguments
%
% subject_lengths -- First sequence of collinear points
%   A vector of length 'm' containing the 1D coordinates of points on a
%   line. The elements of `subject_lengths` should be sorted.
%
% query_lengths -- Second sequence of collinear points
%   A vector of length 'n' containing the 1D coordinates of points on the
%   same or a different line from the points in `subject_lengths`. The
%   points in `subject_lengths` and `query_lengths` can be in different
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
% color_weight -- Weight of colour adjacency scores
%   The weight to apply to the colour portion of the alignment scores.
%   `color_weight` ranges from `0` to `1`, where `0` enforces matching
%   based only on cross ratios, and `1` enforces matching based only on
%   colour adjacency constraints.
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
% The two series of points are matched based on their cross ratios, because
% cross ratios (as opposed to positions or length ratios) are
% invariant to projective transformations:
% 
% First, all possible cross ratios (from all possible combinations of
% points) are generated for the two sets of points. Each cross ratio from
% the subject points is compared to each cross ratio from the query points,
% and the degree of agreement between the two cross ratios is computed. The
% agreement score is then added to the cross ratio matching scores of the
% points used to form the two cross ratios.
%
% For example, if the cross ratios of points `(a, b, c, d)` and `(i, j, k,
% l)` were compared, the resulting score would be added to the matching
% scores for the pairings `(a, i)`, `(b, j)`, `(c, k)`, and `(d, l)`.
%
% The cross ratio-based matching scores are linearly rescaled to the range
% `[0,1]` using the minimum and maximum of the scores.
%
% Second, if colour arguments have been provided, a second table of scores
% is computed to reflect the agreement between the colours to the left and
% right of the points.
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
% The final matching scores for the points in the subject and query
% sequences are weighted averages of the cross ratio-based scores and the
% colour-based scores. The weight on the colour-based scores is
% `color_weight`.
%
% Third, the sequences of points are aligned using
% `swSequenceAlignmentAffine`, in combination with the matching scores
% determined as described above. A semi-global alignment scheme is used, as
% `query` is assumed to contain valid points, but possibly be missing
% points, whereas `subject` is assumed to be the sequence of all valid
% points. Consequently, the complete alignment of `query` will be favoured,
% whereas the process will not penalize the extension of `subject` outside
% the aligned region.
%
% See also crossRatio, swSequenceAlignmentAffine

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created February 9, 2017


end

