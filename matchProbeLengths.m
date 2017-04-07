function [ subject_match_indices ] = matchProbeLengths(...
    subject_lengths, query_lengths,...
    subject_gap_cost, query_gap_cost,...
    affine_weight, varargin...
)
% MATCHPROBELENGTHS  Align sequences of points on two lines using cross ratios and colour adjacency relationships
%
% ## Syntax
% subject_match_indices = matchProbeLengths(...
%   subject_lengths, query_lengths, subject_gap_cost, query_gap_cost,...
%   affine_weight [, verbose]...
% )
%
% subject_match_indices = matchProbeLengths(...
%   subject_lengths, query_lengths, subject_gap_cost, query_gap_cost,...
%   affine_weight, subject_colors, query_colors, color_weight...
%   [, verbose]...
% )
%
% ## Description
% subject_match_indices = matchProbeLengths(...
%   subject_lengths, query_lengths, subject_gap_cost, query_gap_cost,...
%   affine_weight [, verbose]...
% )
%   Returns the indices of points in the subject sequence which match
%   points in the query sequence. The matching is based on point cross
%   ratios and/or length ratios.
%
% subject_match_indices = matchProbeLengths(...
%   subject_lengths, query_lengths, subject_gap_cost, query_gap_cost,...
%   affine_weight, subject_colors, query_colors, color_weight...
%   [, verbose]...
% )
%   Returns the indices of points in the subject sequence which match
%   points in the query sequence. The matching is based on point cross
%   ratios and/or length ratios, and on the colours to the left and right
%   of each point.
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
% affine_weight -- Weight of length ratio scores
%   Under affine distortion conditions, points can be aligned using length
%   ratios, whereas under projective distortion conditions, length ratios
%   are no longer preserved, and points should be aligned using cross
%   ratios.
%
%   An affine approximation is relatively accurate in many cases, and
%   length ratios tend to be more discriminating than cross ratios.
%   Furthermore, length ratios are faster to compute.
%
%   This parameter is the weight to apply to the length ratio-based portion
%   of the alignment scores with respect to the cross ratio-based portion
%   of the alignment scores. It ranges from `0` to `1`, where `0` enforces
%   matching based on cross ratios, and `1` enforces matching based on
%   length ratios.
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
%   based only on geometric constraints (cross ratios and/or length
%   ratios), and `1` enforces matching based only on colour adjacency
%   constraints.
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
% The two series of points are matched based on a combination of their
% cross ratios, because cross ratios (as opposed to positions or length
% ratios) are invariant to projective transformations, and their length
% ratios, because affine distortion is often a reasonable approximation to
% the actual transformation.
% 
% First, all possible cross ratios (from all possible 4-point subsequences)
% are generated for the two sets of points. Each cross ratio from the
% subject points is compared to each cross ratio from the query points, and
% the degree of agreement between the two cross ratios is computed. The
% agreement score is then added to the cross ratio matching scores of the
% points used to form the two cross ratios.
%
% The same procedure is repeated for length ratios, except using all
% possible 3-point subsequences.
%
% For example, if the cross ratios of points `(a, b, c, d)` and `(i, j, k,
% l)` were compared, the resulting score would be added to the matching
% scores for the pairings `(a, i)`, `(b, j)`, `(c, k)`, and `(d, l)`.
%
% The cross ratio-based matching scores are normalized to the range `[0,1]`
% using the counts of cross ratios compared to generate each score. The
% same normalization procedure is applied to length ratio-based matching
% scores, and the two sets of scores are combined, using the
% `affine_weight` weighting factor, to form the final geometric matching
% scores.
%
% Next, if colour arguments have been provided, a second table of scores
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
% sequences are weighted averages of the geometry-based scores and the
% colour-based scores. The weight on the colour-based scores is
% `color_weight`.
%
% Finally, the sequences of points are aligned using
% `swSequenceAlignmentAffine`, in combination with the matching scores
% determined as described above. A semi-global alignment scheme is used, as
% `query` is assumed to contain valid points, but possibly be missing
% points, whereas `subject` is assumed to be the sequence of all valid
% points. Consequently, the complete alignment of `query` will be favoured,
% whereas the process will not penalize the extension of `subject` outside
% the aligned region.
%
% The alignment is tested in both the forward and reverse directions, and
% the final output, `subject_match_indices`, corresponds to the direction
% producing the highest alignment score. Note that geometric and
% colour-based point matching scores must be computed for each of the
% alignment directions.
%
% See also crossRatio, lengthRatio, swSequenceAlignmentAffine

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created February 9, 2017

    function [ score ] = f_similarity_forward(s, q)
        score = scores(s, q, 1);
    end

    function [ score ] = f_similarity_reverse(s, q)
        score = scores(s, q, 2);
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
narginchk(5, 9);

% Parse input arguments
match_colors = false;
verbose = false;
if ~isempty(varargin)
    if length(varargin) == 1
        verbose = varargin{1};
    elseif length(varargin) >= 3
        match_colors = true;
        color_weight = varargin{3};
        match_colors = match_colors && color_weight > 0;
        subject_colors = varargin{1};
        query_colors = varargin{2};
        if length(varargin) == 4
            verbose = varargin{4};
        end
    else
        error('Incorrect number of input arguments')
    end
else
    
end

n_subject = length(subject_lengths);
subject_sequence = 1:n_subject;
n_query = length(query_lengths);
query_sequence = 1:n_query;

% Computing the cross ratio requires 4 points. Disambiguating between
% sequence orientations requires 5 points, because the cross ratio is
% invariant to a reversal of point order.
if n_subject < 4
    warning('Insufficient points given on the subject line to compute cross ratios.')
elseif n_subject < 5
    warning('Insufficient points given on the subject line to determine forwards/backwards orientation from cross ratios.')
end
if n_subject < 3
    warning('Insufficient points given on the subject line to compute length ratios.')
end
if n_query < 4
    warning('Insufficient points given on the query line to compute cross ratios.')
elseif n_query < 5
    warning('Insufficient points given on the query line to determine forwards/backwards orientation from cross ratios.')
end
if n_query < 3
    warning('Insufficient points given on the query line to compute length ratios.')
end
match_cross_ratios = (n_subject >= 5 && n_query >= 5) && (affine_weight < 1);
match_length_ratios = (n_subject >= 3 && n_query >= 3) && (affine_weight > 0);
if match_colors
    match_cross_ratios = match_cross_ratios && (color_weight < 1);
    match_length_ratios = match_length_ratios && (color_weight < 1);
end
if ~match_cross_ratios && ~match_length_ratios && ~match_colors
    error('Insufficient points given to align sequences using geometric constraints.')
end

if ~all(subject_lengths == sort(subject_lengths))
    error('The `subject_lengths` input argument should be sorted in ascending order.');
end
if ~all(query_lengths == sort(query_lengths))
    error('The `query_lengths` input argument should be sorted in ascending order.');
end

flag_index = 1;
for flag = [match_cross_ratios, match_length_ratios]
    if flag
        if flag_index == 1
            subsequence_length = 4;
            ratioFcn = @crossRatio;
            str = 'cross ratio matching scores';
        else
            subsequence_length = 3;
            ratioFcn = @lengthRatio;
            str = 'length ratio matching scores';
        end
        
        % Compute all possible ratios of the subject and query points
        % Note that 'nchoosek' preserves the order of the items being chosen.
        subject_combinations = nchoosek(subject_sequence, subsequence_length);
        n_subject_ratios = size(subject_combinations, 1);
        subject_ratios = zeros(n_subject_ratios, 1);
        for i = 1:n_subject_ratios
            points = subject_lengths(subject_combinations(i, :), 1);
            subject_ratios(i) = ratioFcn(points);
        end

        query_combinations = nchoosek(query_sequence, subsequence_length);
        n_query_ratios = size(query_combinations, 1);
        query_ratios = zeros(n_query_ratios, 1);
        for i = 1:n_query_ratios
            points = query_lengths(query_combinations(i, :), 1);
            query_ratios(i) = ratioFcn(points);
        end

        % Compute ratio scores for the forward and reverse alignments
        % The first layer contains scores for forward alignment; The second
        % contains scores for reverse alignment.
        ratio_scores = zeros(n_subject, n_query, 2);
        ratio_score_counts = zeros(n_subject, n_query, 2);
        point_indices_query = [
            reshape(query_combinations.', [], 1);
            reshape(fliplr(query_combinations).', [], 1)
            ];
        point_indices_direction = [
            ones(subsequence_length * n_query_ratios, 1);
            2 * ones(subsequence_length * n_query_ratios, 1)
            ];
        % `duplicates_buffer` prevents clobbering that would otherwise occur when
        % the same pairing of a subject and a query point corresponds to multiple
        % ratio scores within an iteration of the cross ratio scoring loop
        % below.
        ratio_position_indices = repmat((1:subsequence_length).', 2 * n_query_ratios, 1); % Index of position in cross ratio
        ratio_score_indices = sub2ind(...
            size(ratio_scores),...
            ratio_position_indices,...
            point_indices_query,...
            point_indices_direction...
            );
        n_duplicates_max = nchoosek(n_query - 1, subsequence_length - 1);
        [unique_indices,unique_indices_map,unique_indices_rows] = unique(ratio_score_indices);
        n_unique_indices = length(unique_indices);
        duplicates_buffer = zeros(n_unique_indices, n_duplicates_max);
        unique_indices_columns = zeros(length(ratio_score_indices), 1);
        for i = 1:n_unique_indices
            unique_index_filter = (unique_indices_rows == i);
            unique_indices_columns(unique_index_filter) = (1:sum(unique_index_filter)).';
        end
        duplicates_buffer_indices = sub2ind(...
            size(duplicates_buffer), unique_indices_rows, unique_indices_columns...
            );
        duplicates_buffer(duplicates_buffer_indices) = 1;
        ratio_score_counts_increment = sum(duplicates_buffer, 2);

        for i = 1:n_subject_ratios
            point_indices_subject = repmat(subject_combinations(i, :)', n_query_ratios * 2, 1);
            ratio_subject = repmat(subject_ratios(i), 2 * n_query_ratios, 1);
            if flag_index == 1
                % Cross ratios are invariant to point sequence reversal
                ratio_score = abs(ratio_subject - repmat(query_ratios, 2, 1)) ./...
                    ratio_subject;
            else
                % Length ratios are inverted when points are reversed
                ratio_score = abs(ratio_subject - [query_ratios; query_ratios.^(-1)]) ./...
                    ratio_subject;
            end
            ratio_score = max(1 - ratio_score, 0);
            ratio_score_indices = sub2ind(...
                size(ratio_scores),...
                point_indices_subject,...
                point_indices_query,...
                point_indices_direction...
                );
            ratio_score = repelem(ratio_score, subsequence_length);
            duplicates_buffer(duplicates_buffer_indices) = ratio_score;
            ratio_score_indices = ratio_score_indices(unique_indices_map);
            ratio_scores(ratio_score_indices) =...
                ratio_scores(ratio_score_indices) +...
                sum(duplicates_buffer, 2);
            ratio_score_counts(ratio_score_indices) =...
                ratio_score_counts(ratio_score_indices) +...
                ratio_score_counts_increment;
        end

        % Normalize
        ratio_scores = ratio_scores ./ ratio_score_counts;
        ratio_scores(isnan(ratio_scores)) = 0;

        if verbose
            plotScores(ratio_scores, str);
        end
        
        if flag_index == 1 || (flag_index == 2 && ~match_cross_ratios)
            combined_ratio_scores = ratio_scores;
        else
            combined_ratio_scores = (ratio_scores * affine_weight) +...
                (combined_ratio_scores * (1 - affine_weight));
            
            if verbose
                plotScores(combined_ratio_scores, 'cross ratio and length ratio combined matching scores');
            end
        end
    end
    flag_index = 2;
end

if match_colors
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
    
    if match_cross_ratios || match_length_ratios
        scores = (color_scores * color_weight) +...
            (combined_ratio_scores * (1 - color_weight));
    
        if verbose
            plotScores(scores, 'combined matching scores');
        end
    else
        scores = color_scores;
    end
else
    scores = combined_ratio_scores;
end

% Match sequences in the given directions with dynamic programming
threshold = -Inf;
subject_gap_cost = [subject_gap_cost 0];
query_gap_cost = [query_gap_cost 0];
[ alignment_forward, score_forward ] = swSequenceAlignmentAffine(...
        subject_sequence, query_sequence,...
        @f_similarity_forward, threshold,...
        subject_gap_cost, query_gap_cost, 'SemiGlobal'...
    );

if verbose
    disp('Forward sequence alignment score:');
    disp(score_forward);
end

% Match sequences in the reverse directions with dynamic programming
query_sequence_reverse = fliplr(query_sequence);
[ alignment_reverse, score_reverse ] = swSequenceAlignmentAffine(...
        subject_sequence, query_sequence_reverse,...
        @f_similarity_reverse, threshold,...
        subject_gap_cost, query_gap_cost, 'SemiGlobal'...
    );

if verbose
    disp('Reverse sequence alignment score:');
    disp(score_reverse);
end

if score_forward > score_reverse
    alignment = alignment_forward;
elseif score_forward == score_reverse
    warning('Scores for forward and reverse alignments are equal. Orientation is ambiguous. Forward orientation is assumed.')
    alignment = alignment_forward;
else
    alignment = alignment_reverse;
    alignment_filter = (alignment(:, 2) ~= 0);
    alignment(alignment_filter, 2) = query_sequence_reverse(alignment(alignment_filter, 2));
end

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

