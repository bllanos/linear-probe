%% Test cases for swSequenceAlignmentAffine

% Bernard Llanos
% Spring 2016 research assistantship supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created August 12, 2016

%% Local sequence alignment
% Example from https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm

% A = 1, C = 2, G = 3, T = 4
subject = [1 3 2 1 2 1 2 1];
query = [1 2 1 2 1 2 4 1];
f_similarity = @(x, y) (x == y) * 2 + (x ~= y) * -1;
threshold = 0;
subject_gap_penalty = [-1 0];
query_gap_penalty = [-1 0];
[ alignment, score ] = swSequenceAlignmentAffine( subject, query, f_similarity, threshold, subject_gap_penalty, query_gap_penalty, 'Local' );
alignment_expected = [
        1     1;
        2     0;
        3     2;
        4     3;
        5     4;
        6     5;
        7     6;
        0     7;
        8     8
    ];
score_expected = 12;
if all(all(alignment_expected == alignment))
    disp('Local alignment 1 is correct.')
else
    disp('Incorrect local alignment 1.')
end
if score_expected == score
    disp('Local alignment 1 score is correct.')
else
    disp('Incorrect local alignment 1 score.')
end

%% Global sequence alignment
% Example from http://pages.cs.wisc.edu/~bsettles/ibs08/lectures/02-alignment.pdf
% (Page 27)

% A = 1, C = 2, G = 3, T = 4
subject = [1 1 1 2];
query = [1 3 2];
f_similarity = @(x, y) (x == y) * 1 + (x ~= y) * -1;
threshold = -Inf;
subject_gap_penalty = [-2 -2];
query_gap_penalty = [-2 -2];
[ alignment, score ] = swSequenceAlignmentAffine( subject, query, f_similarity, threshold, subject_gap_penalty, query_gap_penalty, 'Global' );
alignment_expected = [
        2 1;
        3 2;
        4 3
    ];
score_expected = -1;
if all(all(alignment_expected == alignment))
    disp('Global alignment 1 is correct.')
else
    disp('Incorrect global alignment 1.')
end
if score_expected == score
    disp('Global alignment 1 score is correct.')
else
    disp('Incorrect global alignment 1 score.')
end

%% Semi-global sequence alignment
% Example from http://ocw.mit.edu/courses/electrical-engineering-and-computer-science/6-096-algorithms-for-computational-biology-spring-2005/lecture-notes/lecture5_newest.pdf
% (Page 24)

% A = 1, C = 2, G = 3, T = 4
subject = [1 1 3 2];
query = [1 3 4];
f_similarity = @(x, y) (x == y) * 1 + (x ~= y) * -1;
threshold = -Inf;
subject_gap_penalty = [-2 -2];
query_gap_penalty = [-2 -2];
[ alignment, score ] = swSequenceAlignmentAffine( subject, query, f_similarity, threshold, subject_gap_penalty, query_gap_penalty, 'SemiGlobal' );
alignment_expected = [
        2 1;
        3 2;
        4 3
    ];
score_expected = 1;
if all(all(alignment_expected == alignment))
    disp('Semi-global alignment 1 is correct.')
else
    disp('Incorrect semi-global alignment 1.')
end
if score_expected == score
    disp('Semi-global alignment 1 score is correct.')
else
    disp('Incorrect semi-global alignment 1 score.')
end

%% Semi-global sequence alignment 2
% Example from http://ocw.mit.edu/courses/electrical-engineering-and-computer-science/6-096-algorithms-for-computational-biology-spring-2005/lecture-notes/lecture5_newest.pdf
% (Page 22)

% A = 1, C = 2, G = 3, T = 4
subject = [2 1 3 2 1 2 4 4 3 3 1 4 4 2 4 2 3 3];
query = [2 1 3 2 3 4 3 3];
f_similarity = @(x, y) (x == y) * 1 + (x ~= y) * -1;
threshold = -Inf;
subject_gap_penalty = [-2 -2];
query_gap_penalty = [-2 -2];
[ alignment, score ] = swSequenceAlignmentAffine( subject, query, f_similarity, threshold, subject_gap_penalty, query_gap_penalty, 'SemiGlobal' );
alignment_expected = [
        4 1;
        5 2;
        0 3;
        6 4;
        7 5;
        8 6;
        9 7;
        10 8
    ];
score_expected = 3;
if all(all(alignment_expected == alignment))
    disp('Semi-global alignment 2 is correct.')
else
    disp('Incorrect semi-global alignment 2.')
end
if score_expected == score
    disp('Semi-global alignment 2 score is correct.')
else
    disp('Incorrect semi-global alignment 2 score.')
end

%% Local sequence alignment 2
% Example from http://pages.cs.wisc.edu/~bsettles/ibs08/lectures/02-alignment.pdf
% (Page 35)

% A = 1, C = 2, G = 3, T = 4
subject = [4 4 1 1 3];
query = [1 1 3 1];
f_similarity = @(x, y) (x == y) * 1 + (x ~= y) * -1;
threshold = 0;
subject_gap_penalty = [-2 -2];
query_gap_penalty = [-2 -2];
[ alignment, score ] = swSequenceAlignmentAffine( subject, query, f_similarity, threshold, subject_gap_penalty, query_gap_penalty, 'Local' );
alignment_expected = [
        3 1;
        4 2;
        5 3
    ];
score_expected = 3;
if all(all(alignment_expected == alignment))
    disp('Local alignment 2 is correct.')
else
    disp('Incorrect local alignment 2.')
end
if score_expected == score
    disp('Local alignment 2 score is correct.')
else
    disp('Incorrect local alignment 2 score.')
end

%% Global sequence alignment with affine gap function
% Example from http://pages.cs.wisc.edu/~bsettles/ibs08/lectures/02-alignment.pdf
% (Page 41)

% A = 1, C = 2, G = 3, T = 4
subject = [1 1 4];
query = [1 2 1 2 4];
f_similarity = @(x, y) (x == y) * 1 + (x ~= y) * -1;
threshold = -Inf;
subject_gap_penalty = [-4 -1];
query_gap_penalty = [-4 -1];
[ alignment, score ] = swSequenceAlignmentAffine( subject, query, f_similarity, threshold, subject_gap_penalty, query_gap_penalty, 'Global' );
alignment_expected = [
        1 3;
        2 4;
        3 5
    ];
score_expected = -4;
if all(all(alignment_expected == alignment))
    disp('Global affine gap penalty alignment 1 is correct.')
else
    disp('Incorrect global affine gap penalty alignment 1.')
end
if score_expected == score
    disp('Global affine gap penalty alignment 1 score is correct.')
else
    disp('Incorrect global affine gap penalty alignment 1 score.')
end