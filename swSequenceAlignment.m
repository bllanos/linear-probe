function [ alignment, score ] = swSequenceAlignment( subject, query, f_similarity, threshold, f_subject_gap, f_query_gap, type )
% SWSEQUENCEALIGNMENT  Optimal global, semi-global, or local pairwise sequence alignment
%
% ## Syntax
% alignment = swSequenceAlignment( subject, query, f_similarity, threshold, f_subject_gap, f_query_gap, type )
% [ alignment, score ] = swSequenceAlignment( subject, query, f_similarity, threshold, f_subject_gap, f_query_gap, type )
%
% ## Description
% alignment = swSequenceAlignment( subject, query, f_similarity, threshold, f_subject_gap, f_query_gap, type )
%   Returns the overlapping region of the subject and query sequences in an
%   optimal alignment.
%
% [ alignment, score ] = swSequenceAlignment( subject, query, f_similarity, threshold, f_subject_gap, f_query_gap, type )
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
% f_subject_gap -- Score function for gaps in the subject sequence
%   The handle to a function of one integer value (class 'double') which
%   returns a numeric scalar output argument. `f_subject_gap` is passed the
%   number of consecutive positions for which an element of the query
%   sequence is matched with a gap in the subject sequence. It should
%   return the score corresponding to the gap length.
%
% f_query_gap -- Score function for gaps in the query sequence
%   Similar to `f_subject_gap`, but returns scores for situations in which
%   elements in the subject sequence are matched with gaps in the query
%   sequence.
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
%     using `f_subject_gap` and `f_query_gap`. If gaps are penalized with
%     low scores, then this encourages the alignment to span the length of
%     both sequences.
%   - 'SemiGlobal', then leading and trailing gaps in the query sequence
%     are scored using `f_query_gap`, but leading and trailing gaps in the
%     subject sequence are given scores of zero, independent of their
%     lengths. If gaps are penalized with low scores, then this encourages
%     the alignment to span the length of the query sequence, but not
%     necessarily the subject sequence.
%   - 'Local', then leading and trailing gaps in either sequence are given
%     scores of zero, independent of their lengths. If gaps are penalized
%     with low scores, then this scheme makes it more likely that the final
%     alignment will span only a portion of each sequence.
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
% score -- Optimal sequence alignment score
%   The score corresponding to the alignment returned in `alignment`, which
%   is the optimal (maximum) alignment score.
%
% ## Notes
% - This is, essentially, a naive implementation of the dynamic programming
%   Smith-Waterman algorithm (see references below).
% - The gap score functions `f_subject_gap` and `f_query_gap` can be
%   general convex gap penalty functions, including affine gap penalty
%   functions. The only limitation is that the gap scores depend only on
%   the length of the gap and not on the position of the gap or other
%   contextual information.
%
% ## References
% - Primary reference. If `threshold` is zero, and `type` is `Local`, then
%   the algorithm behaves as presented at:
%   https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm
% - Reference used for algorithm behaviour in the global and semi-global
%   alignment cases:
%   http://ocw.mit.edu/courses/electrical-engineering-and-computer-science/6-096-algorithms-for-computational-biology-spring-2005/lecture-notes/lecture5_newest.pdf

% Bernard Llanos
% Spring 2016 research assistantship supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created August 3, 2016

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

% Initialization
m = length(subject);
n = length(query);
H = zeros(m + 1, n + 1);
back_pointers = zeros(m + 1, n + 1, 2);

% Compute possible gap penalties
subject_gaps = arrayfun(f_subject_gap, 1:n);
query_gaps = arrayfun(f_query_gap, 1:m);

if strcmp(type, 'Global')
    H(1, 2:end) = subject_gaps;
    back_pointers(1, 2:end, 1) = 1;
    back_pointers(1, 2:end, 2) = 1:n;
    H(2:end, 1) = query_gaps;
    back_pointers(2:end, 1, 1) = 1:m;
    back_pointers(2:end, 1, 2) = 1;
elseif strcmp(type, 'SemiGlobal')
    H(1, 2:end) = subject_gaps;
    back_pointers(1, 2:end, 1) = 1;
    back_pointers(1, 2:end, 2) = 1:n;
end

% Compute best alignment score
for i = 2:(m + 1)
    query_gaps_i = query_gaps(1:(i-1));
    for j = 2:(n + 1)
        subject_gaps_j = subject_gaps(1:(j-1));
        [value_gap_subject, length_gap_subject] = max(H(i, (j-1):-1:1) + subject_gaps_j);
        [value_gap_query, length_gap_query] = max(H((i-1):-1:1, j).' + query_gaps_i);
        value_match = H(i-1, j-1) + f_similarity(subject(i-1), query(j-1));
        [H(i, j), ind] = max([threshold, value_match, value_gap_query, value_gap_subject]);
        switch ind
            case 2
                back_pointers(i, j, :) = [i-1, j-1];
            case 3
                back_pointers(i, j, :) = [i-length_gap_query, j];
            case 4
                back_pointers(i, j, :) = [i, j-length_gap_subject];
        end
    end
end

% Trace backpointers to recover the path
path = zeros(m + n, 2);
if strcmp(type, 'Global')
    score = H(end);
    i = m + 1;
    j = n + 1;
elseif strcmp(type, 'SemiGlobal')
    [score, i] = max(H(:, end));
    j = n + 1;
else
    % Local alignment
    [column_scores, row_indices] = max(H);
    [score, j] = max(column_scores);
    i = row_indices(j);
end

path(1, :) = [i, j];
k = 1;
while i > 1 && j > 1
    k = k + 1;
    path(k, :) = back_pointers(i, j, :);
    i = path(k, 1);
    j = path(k, 2);
end

% Trace the path to recover the alignment
path = [ ones(1, 2); path(k:-1:1, :) ];
path_diff = diff(path);
alignment = zeros(k, 2);
for s = 1:k
    pair = path(s + 1, :) - 1;
    if path_diff(s, 1)
        if path_diff(s, 2)
            % Match or mismatch
            alignment(s, :) = pair;
        else
            % Gap in query
            alignment(s, 1) = pair(1);
        end
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