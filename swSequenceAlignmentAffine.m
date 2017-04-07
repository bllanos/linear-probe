function [ alignment, score ] = swSequenceAlignmentAffine(...
    subject, query, f_similarity, threshold,...
    subject_gap_penalty, query_gap_penalty, type, varargin...
)
% SWSEQUENCEALIGNMENT  Optimal global, semi-global, or local pairwise sequence alignment for affine or simpler gap penalty functions
%
% ## Syntax
% alignment = swSequenceAlignmentAffine(...
%   subject, query, f_similarity, threshold,...
%   subject_gap_penalty, query_gap_penalty, type [, sample_count]...
% )
% [ alignment, score ] = swSequenceAlignmentAffine(...
%   subject, query, f_similarity, threshold,...
%   subject_gap_penalty, query_gap_penalty, type [, sample_count]...
% )
%
% ## Description
% alignment = swSequenceAlignmentAffine(...
%   subject, query, f_similarity, threshold, subject_gap_penalty, query_gap_penalty, type...
% )
%   Returns the overlapping region of the subject and query sequences in an
%   optimal alignment.
%
% [ alignment, score ] = swSequenceAlignmentAffine(...
%   subject, query, f_similarity, threshold, subject_gap_penalty, query_gap_penalty, type...
% )
%   Additionally returns the score of the optimal alignment.
%
% ## Input Arguments
%
% subject -- First sequence
%   A numeric vector of length 'm' containing the first sequence in the
%   pair. The "subject" sequence is the sequence to be partially aligned,
%   in the case of semi-global alignment. (In a semi-global alignment, the
%   overlap region need not span the length of the subject sequence.)
%
% query -- Second sequence
%   A numeric vector of length 'n' containing the second sequence in the
%   pair. The "query" sequence is the sequence to be fully-aligned,
%   in the case of semi-global alignment. (In a semi-global alignment, the
%   overlap region will span the length of the query sequence.)
%
% f_similarity -- Sequence element match/mismatch score function
%   The handle to a function accepting two scalar numeric input arguments
%   (of the same classes as 'subject' and 'query', respectively) and
%   returning a numeric scalar output argument. `f_similarity(a, b)`
%   returns the score of a match between an element `a` from `subject` and
%   an element `b` from `query`. A higher score indicates a more suitable
%   match.
%
% threshold -- Local alignment termination value
%   The alignment score below which the algorithm will stop extending the
%   current alignment and start a new local alignment. `threshold` allows
%   for local sequence alignment in which poorly-matching regions of the
%   subject and query sequences are excluded. It is to be used in
%   conjunction with a `type` input argument of 'Local'.
%
%   The appropriate value of `threshold` depends on the scores given to
%   gaps and element matches. For global or semi-global alignments,
%   `threshold` should be `-Inf`.
%
% subject_gap_penalty -- Scores for gaps in the subject sequence
%   A two-element vector containing the scores for situations in which
%   elements in the query sequence are matched with gaps in the subject
%   sequence. The first element is the score for a gap of length 1. The
%   second element is the score for each additional increment in the gap
%   length.
%
%   This construction allows for constant, linear, or affine gap penalty
%   functions.
%
%   Note that the first element of `subject_gap_penalty` is not the cost
%   for starting a gap (which is the cost of a gap of length zero), in
%   contrast to the convention used in the references. The cost of a gap of
%   length zero is calculated by this function as
%   `subject_gap_penalty(1) - subject_gap_penalty(2)`.
%
% query_gap_penalty -- Scores for gaps in the query sequence
%   Similar to `subject_gap_penalty`, but contains scores for situations in
%   which elements in the subject sequence are matched with gaps in the
%   query sequence.
%
% type -- Alignment scheme
%   A string equal to 'Global', 'SemiGlobal', or 'Local', which describes
%   the treatment of leading and trailing gaps. A leading gap is a
%   situation in which the start of one sequence is aligned with gaps,
%   whereas a trailing gap is a situation in which the end of one sequence
%   is aligned with gaps.
%
%   If `type` is:
%   - 'Global', then leading and trailing gaps in either sequence are scored
%     using `subject_gap_penalty` and `query_gap_penalty`. If gaps are
%     penalized with low scores, then this encourages the alignment to span
%     the length of both sequences.
%   - 'SemiGlobal', then leading and trailing gaps in the query sequence
%     are scored using `query_gap_penalty`, but leading and trailing gaps
%     in the subject sequence are given scores of zero, independent of
%     their lengths. If gaps are penalized with low scores, then this
%     encourages the alignment to span the length of the query sequence,
%     but not necessarily the subject sequence.
%   - 'Local', then leading and trailing gaps in either sequence are given
%     scores of zero, independent of their lengths. If gaps are penalized
%     with low scores, then this scheme makes it more likely that the final
%     alignment will span only a portion of each sequence.
%
% sample_count -- Number of alignments to output
%   If `sample_count` is greater than one, `alignment` will be a cell row
%   vector of length `sample_count`. Each element of `alignment` will be an
%   optimal sequence alignment sampled using a random walk through the
%   forest of back pointers defining possible optimal sequence alignments.
%
%   `sample_count` is used to obtain a representative set of optimal
%   alignments in cases where the optimal alignment is not unique.
%   Enumerating every possible optimal alignment may be computationally
%   infeasible; A fixed number of samples trades completeness for
%   efficiency.
%
%   Defaults to 1 if not passed.
%
% ## Output Arguments
%
% alignment -- Optimal sequence alignment
%   The optimal global, semi-global, or local sequence alignment (depending
%   on the value of `type`) obtained using the Smith-Waterman algorithm.
%   The alignment is represented as a k x 2 array, where `alignment(i, :)`
%   indicates that `subject(alignment(i, 1))` is matched with
%   `query(alignment(i, 2))`. If `alignment(i, j)` is zero, it indicates
%   that, if `j == 1`, `query(alignment(i, 2))` is matched with a gap, or,
%   if `j == 2`, `subject(alignment(i, 2))` is matched with a gap.
%   
%   Leading and trailing gaps are not output in `alignment`, as they can be
%   inferred. For example:
%   - If `alignment(1, :)` is `[ 3 1 ]`, then `subject(1:2)` is matched
%     with a leading gap in `query`.
%   - If `alignment(end, :)` is `[ 15 12 ]`, then either `subject` has
%     length 15, and `query(13:end)` is matched with a trailing gap in
%     `subject`, or `query` has length 12, and `subject(15:end)` is matched
%     with a trailing gap in `query`.
%
%   If `sample_count` is greater than 1, `alignment` is a cell vector of
%   optimal sequence alignments. Refer to the documentation of
%   `sample_count` above.
%
% score -- Optimal sequence alignment score
%   The score corresponding to the alignment returned in `alignment`, which
%   is the optimal (maximum) alignment score.
%
% ## Notes
% - This is an implementation of the dynamic programming
%   Smith-Waterman algorithm for an affine gap penalty (see references
%   below).
%
% ## References
% - Primary reference. If `threshold` is zero, and `type` is `Local`, then
%   the algorithm behaves as presented at:
%   https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm
% - Reference used for algorithm behaviour in the global and semi-global
%   alignment cases, and for affine gap penalty functions:
%   http://ocw.mit.edu/courses/electrical-engineering-and-computer-science/6-096-algorithms-for-computational-biology-spring-2005/lecture-notes/lecture5_newest.pdf
% - Additional reference used for algorithm behaviour with respect to an
%   affine gap penalty function:
%   http://pages.cs.wisc.edu/~bsettles/ibs08/lectures/02-alignment.pdf
%
% See also swSequenceAlignmentAffine

% Bernard Llanos
% Spring 2016 research assistantship supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created August 12, 2016

nargoutchk(1, 2);
narginchk(7, 7);

if ~(...
        strcmp(type, 'Global') ||...
        strcmp(type, 'SemiGlobal') ||...
        strcmp(type, 'Local')...
    )
    error('Unknown value of `type` input argument');
end

if ~strcmp(type, 'Local') && threshold ~= -Inf
    warning('For a global or semi-global alignment, the value of `threshold` should be `-Inf`.');
end

if ~isempty(varargin)
    sample_count = varargin{1};
else
    sample_count = 1;
end
if sample_count < 1
    error('Value of `sample_count` input argument must be one or greater.')
end

% Initialization
m = length(subject);
n = length(query);
% In http://pages.cs.wisc.edu/~bsettles/ibs08/lectures/02-alignment.pdf,
% page 38:
% - H(:, :, 1) is 'M'
% - H(:, :, 2) is 'I_x'
% - H(:, :, 3) is 'I_y'
H = -Inf(m + 1, n + 1, 3);
back_pointers = zeros(m + 1, n + 1, 3, 3);
subject_gap_penalty_1 = subject_gap_penalty(1);
subject_gap_penalty_2 = subject_gap_penalty(2);
query_gap_penalty_1 = query_gap_penalty(1);
query_gap_penalty_2 = query_gap_penalty(2);

if strcmp(type, 'Global')
    H(1, :, 3) = subject_gap_penalty_1 + subject_gap_penalty_2 * (-1:1:n-1);
    back_pointers(1, 2:end, 3, 1) = 1;
    back_pointers(1, 2:end, 3, 2) = 1:n;
    back_pointers(1, 2:end, 3, 3) = 3;
    H(:, 1, 2) = query_gap_penalty_1 + query_gap_penalty_2 * (-1:1:m-1);
    back_pointers(2:end, 1, 2, 1) = 1:m;
    back_pointers(2:end, 1, 2, 2) = 1;
    back_pointers(2:end, 1, 2, 3) = 2;
    H(1, 1, 1) = 0;
elseif strcmp(type, 'SemiGlobal')
    H(1, :, 3) = subject_gap_penalty_1 + subject_gap_penalty_2 * (-1:1:n-1);
    back_pointers(1, 2:end, 3, 1) = 1;
    back_pointers(1, 2:end, 3, 2) = 1:n;
    back_pointers(1, 2:end, 3, 3) = 3;
    H(:, 1, 1) = zeros(m + 1, 1);
else
    % Local alignment
    H(1, :, 1) = zeros(1, n + 1);
    H(:, 1, 1) = zeros(m + 1, 1);
end

% Compute best alignment score
for i = 2:(m + 1)
    for j = 2:(n + 1)
        
        % Update 'M'
        value_match = f_similarity(subject(i-1), query(j-1));
        value_match_previous = H(i-1, j-1, 1) + value_match;
        value_gap_query = H(i-1, j-1, 2) + value_match;
        value_gap_subject = H(i-1, j-1, 3) + value_match;
        [H(i, j, 1), ind] = max([value_match_previous, value_gap_query, value_gap_subject, threshold]);
        if ind ~= 4
            back_pointers(i, j, 1, :) = [i-1, j-1, ind];
        end
        
        % Update 'I_x'
        value_match_previous = H(i-1, j, 1) + query_gap_penalty_1;
        value_gap_query = H(i-1, j, 2) + query_gap_penalty_2;
        value_gap_subject_previous = H(i-1, j, 3) + query_gap_penalty_1;
        [H(i, j, 2), ind] = max([value_match_previous, value_gap_query, value_gap_subject_previous]);
        back_pointers(i, j, 2, :) = [i-1, j, ind];
        
        % Update 'I_y'
        value_match_previous = H(i, j-1, 1) + subject_gap_penalty_1;
        value_gap_query_previous = H(i, j-1, 2) + subject_gap_penalty_1;
        value_gap_subject = H(i, j-1, 3) + subject_gap_penalty_2;
        [H(i, j, 3), ind] = max([value_match_previous, value_gap_query_previous, value_gap_subject]);
        back_pointers(i, j, 3, :) = [i, j-1, ind];
    end
end

    function [alignment, score] = tracePath()
        % Trace backpointers to recover the path
        path = zeros(m + n, 3);
        if strcmp(type, 'Global')
            [score, k] = max(H(end, end, :));
            p = m + 1;
            q = n + 1;
        elseif strcmp(type, 'SemiGlobal')
            [column_scores, row_indices] = max(H(:, end, :));
            [score, k] = max(column_scores);
            p = row_indices(k);
            q = n + 1;
        else
            % Local alignment
            [sheet_column_scores, sheet_row_indices] = max(H);
            [sheet_scores, sheet_column_indices] = max(sheet_column_scores);
            [score, k] = max(sheet_scores);
            q = sheet_column_indices(k);
            p = sheet_row_indices(1, q, k);
        end

        path(1, :) = [p, q, k];
        path_length = 1;
        while p > 1 && q > 1 && k > 0
            path_length = path_length + 1;
            path(path_length, :) = back_pointers(p, q, k, :);
            p = path(path_length, 1);
            q = path(path_length, 2);
            k = path(path_length, 3);
        end

        % Trace the path to recover the alignment
        path = [ ones(1, 3); path(path_length:-1:1, :) ];
        alignment = zeros(path_length, 2);
        for s = 1:path_length
            pair = path(s + 1, 1:2) - 1;
            k = path(s + 1, 3);
            if k == 1
                % Match or mismatch
                alignment(s, :) = pair;
            elseif k == 2
                % Gap in query
                alignment(s, 1) = pair(1);
            else
                % Gap in subject
                alignment(s, 2) = pair(2);
            end
        end

        % Strip redundant gap
        if alignment(1, 1) == 0 || alignment(1, 2) == 0
            alignment = alignment(2:end, :);
        end
    end

[alignment, score] = tracePath();

end