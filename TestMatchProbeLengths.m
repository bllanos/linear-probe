%% Test cases for matchProbeLengths

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created February 9, 2017

%% Correct match in forward direction
subject_lengths = [0 1 2 5 8]';
query_lengths = [0 1 2 5 8]';
subject_colors = [
    0 1;
    1 2;
    2 1;
    1 2;
    2 0
    ];
query_colors = [
    0 1;
    1 2;
    2 1;
    1 2;
    2 0
    ];

subject_gap_cost_detection = -0.1;
query_gap_cost_detection = 0.5;
color_weight_detection = 0.5;
verbose = true;

matches = matchProbeLengths(...
  subject_lengths,...
  query_lengths,...
  subject_gap_cost_detection,...
  query_gap_cost_detection,...
  subject_colors,...
  query_colors,...
  color_weight_detection,...
  verbose...
)

%% Correct match in reverse direction
subject_lengths = [0 1 2 5 8]';
query_lengths = [0 3 6 7 8]';
subject_colors = [
    0 1;
    1 2;
    2 1;
    1 2;
    2 0
    ];
query_colors = [
    0 2;
    2 1;
    1 2;
    2 1;
    1 0
    ];

subject_gap_cost_detection = -0.1;
query_gap_cost_detection = 0.5;
color_weight_detection = 0.5;
verbose = true;

matches = matchProbeLengths(...
  subject_lengths,...
  query_lengths,...
  subject_gap_cost_detection,...
  query_gap_cost_detection,...
  subject_colors,...
  query_colors,...
  color_weight_detection,...
  verbose...
)