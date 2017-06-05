%% Probe Object Detection Model Creation
% Associate user-provided measurements of the probe with points marked on
% the image of the probe, such that the image colours can be used to detect
% the probe in other images.
%
% ## Usage
%   Modify parameters and paths to input data in the first code section
%   below, then run.
%
% ## Probe Design Assumptions
%   The probe is a cylindrically-symmetric object composed of conical or
%   cylindrical segments. Ideally, it is perfectly straight. (Some thin
%   objects, such as wooden rods, may appear straight, but actually be
%   curved and therefore not cylindrically symmetrical.) It may taper to a
%   point at one or both ends, although only one end will be used as the
%   active tip.
%
%   The length of the probe should be divided into different-coloured
%   bands, and the pattern of bands should be asymmetrical, such that it is
%   possible to uniquely determine the orientation of the probe, even if
%   some of the probe is occluded. Specifically, asymmetry of band lengths
%   is assumed, as opposed to colour pattern asymmetry. Band edges should
%   be perpendicular to the probe's axis of cylindrical symmetry. It should
%   not be possible for the probe to self-occlude any of the edges between
%   bands.
%
%   At each junction between two coloured bands, the coloured bands should
%   have the same widths.
%
%   Note: All bands should have strong colours. Those with unsaturated
%   colours can be given colour labels of zero, such that they will be
%   ignored. Unsaturated colours cannot be reliably detected, and will also
%   interfere with the correct detection of the other colours.
%
% ## Input
%
% ### Probe measurements
% A '.mat' file containing a structure called 'probe' with the following
% fields:
% - lengths: A vector of distances of edges of bands from the active end of
%     the probe, including a distance of 0.0 for the active end of the
%     probe, and a distance for the other end of the probe (i.e. the length
%     of the entire probe). Distances are measured along the probe's axis
%     of cylindrical symmetry, and are listed in order starting from the
%     active end of the probe.
% - colors: A vector with `length(lengths) - 1` elements specifying colour
%     indices for the bands of the probe. Colour indices for colours that
%     are to be detected should be consecutive integers starting at 1.
%     Colours given indices of zero will be ignored. Indices indicate how
%     the bands of the probe are grouped based on mutually-distinguishable
%     colours, allowing bands to have non-unique colours. The specific
%     index assigned to a given band is unimportant. For instance, the
%     first band need not have a colour index of 1. The elements of
%     `colors` should correspond to adjacent pairs of elements in `lengths`
%     (i.e. The colour indices should be in order starting from the active
%     end of the probe).
% - widths: Width (diameter) of the probe at the edges of bands. Widths
%     must include the ends of the probe, with values of zero or
%     approximate tip diameters for ends that taper to points. (Probe tip
%     widths are currently unused, so their values do not matter.) The
%     elements of 'widths' should correspond to elements of 'lengths'.
%
% Units are arbitrary, but should be consistent with the units of any
% partial reconstruction of an object that the probe is used to refine, and
% with the units used for camera calibration.
%
% ### Image of probe
% An RGB image in any format that can be loaded by `imread` showing the
% entire portion of the probe that the user has provided measurements for.
% The image should depict the probe under lighting conditions that are
% similar to the detection scenario, to facilitate colour-based detection.
%
% The image should not contain significant distortion, or it will be
% difficult to automatically determine the orientation (forwards vs.
% backwards) of the probe in the image.
%
% ### Annotations for probe image
% An image in any format that can be loaded by `imread` with an alpha
% channel that is nonzero at user-marked interest points, and zero
% everywhere else. The image should be based on the same image of the probe
% provided above.
%
% Interest points can be marked by nonzero alpha channel regions that are
% more than one pixel in size - A morphological shrink operation will be
% used to extract single pixel locations from them. Optionally (depending
% on the parameter variable `annotation_corner_search_width` below), a
% corner feature detector will then be used to refine the locations of
% interest points within a search window.
%
% A single interest point should be marked for each end of the probe that
% tapers to a point. Two interest points should be marked for the edges of
% coloured bands on the probe, corresponding to their intersections with
% the probe's contour in the image.
%
% ### Colour noise parameters
% A '.mat' file containing a 'rgb_sigma_polyfit' variable, as output by the
% script '.\EstimateRGBStandardDeviations.m'. 'rgb_sigma_polyfit' describes
% the variation in RGB channel standard deviations with RGB values in the
% image. This information should be computed from images taken under the
% same conditions and with the same camera parameters as the image used for
% probe modeling, if not computed from the same image used for probe
% modeling.
%
% ## Output
%
% ### Probe detection model
% A '.mat' file containing the following variables:
% - 'probe': A copy of the 'probe' variable loaded from the probe
%   measurements file, for convenience.
% - 'probe_band_color_distributions': Discretized variable kernel density
%   estimators of image hue values corresponding to the different coloured
%   bands on the probe, in the same order (starting from the active tip of
%   the probe). The i-th column of this 2D array stores the estimator for
%   the i-th colour class of probe segments.
% - 'probe_band_color_distribution_increment': A scalar equal to the spacing
%   between the samples of hue values in the range [0, 1] at which the
%   variable kernel density estimators have been evaluated to produce
%   'probe_band_color_distributions'. To find the approximate value of the
%   i-th estimator at a query value 'x' in the range [0, 1], use:
%
%   ```
%     queryDiscretized1DFunction(...
%       x, probe_band_color_distributions(:, i),...
%       probe_band_color_distribution_increment...
%     )
%   ```
%
% Additionally, the output file contains the values of all parameters in
% the first section of the script below, for reference. (Specifically,
% those listed in `parameters_list`, which should be updated if the set of
% parameters is changed.)
%
% ## References
% - T. Gevers and H. Stokman. "Robust Histogram Construction from Color
%   Invariants for Object Recognition". IEEE Transactions on Pattern
%   Analysis and Machine Intelligence, vol. 26, no. 1, pp. 113-118, Jan.
%   2004.

% Bernard Llanos
% Spring 2016 research assistantship supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created July 26, 2016

%% Input data and parameters

% List of parameters to save with results
parameters_list = {
        'model_filename',...
        'I_filename',...
        'I_annotations_filename',...
        'annotation_corner_search_width',...
        'point_alignment_outlier_threshold',...
        'subject_gap_cost',...
        'query_gap_cost',...
        'affine_weight',...
        'probe_color_distribution_resolution'
    };

% Probe measurements
model_filename = 'C:\Users\llanos\Google Drive\PointProbing\Data and results\20170410_redPenWithTape\redPenMeasurements.mat';
% Image of probe
I_filename = 'C:\Users\llanos\Google Drive\PointProbing\Data and results\20170410_redPenWithTape\redPenModel.bmp';
% Annotations for image of probe
I_annotations_filename = 'C:\Users\llanos\Google Drive\PointProbing\Data and results\20170410_redPenWithTape\redPenModel_annotated.png';

% Annotation extraction parameters
annotation_corner_search_width = 0; % Set to zero to use centers of user-marked annotations as opposed to nearby corner features

% Parameters for interpreting annotated points
point_alignment_outlier_threshold = 5;

% Parameters for matching annotated points with the probe measurements
subject_gap_cost = -0.1;
query_gap_cost = 0;
affine_weight = 0;

% Number of points at which to evaluate colour estimators
probe_color_distribution_resolution = 180;

% Debugging tools
display_original_image = false;
display_annotations_image = false;
display_extracted_annotations = false;
display_model_from_image = false;
verbose_point_sequence_matching = false;
display_probe_band_masks = false;
display_probe_color_masks = false;
display_hue_image = false;
plot_hue_estimators = true;
plot_hue_classifier = true;

%% Load images and obtain adjusted centers of user-marked annotations

I = imread(I_filename);
image_width = size(I, 2);
image_height = size(I, 1);
image_n_channels = size(I, 3);
if image_n_channels ~= 3
    error('A 3-channel RGB image is required.')
end
if display_original_image
    figure; %#ok<UNRCH>
    imshow(I);
    title('Base image');
end

[~, ~, Ialpha] = imread(I_annotations_filename);
if display_annotations_image
    figure; %#ok<UNRCH>
    imshow(Ialpha);
    title('User-marked annotations');
end

if annotation_corner_search_width
    [ interest_points, feature_search_windows, corners ] = extractInterestPoints( Ialpha, I, annotation_corner_search_width );
else
    interest_points = extractInterestPoints( Ialpha ); %#ok<UNRCH>
end
if display_extracted_annotations
    figure; %#ok<UNRCH>
    I_grey = rgb2gray(I);
    imshow(max(I_grey, Ialpha));
    hold on
    plot(interest_points(:, 1), interest_points(:, 2), 'ro')
    if annotation_corner_search_width
        title('Final marked interest points (red o) and other corner features');
        for rect = feature_search_windows'
            rectangle('Position', rect, 'EdgeColor', 'b');
        end
        for i = 1:length(corners)
            if ~isempty(corners{i})
                plot(corners{i});
            end
        end
    else
        title('Extracted interest points (red o)');
    end
    hold off
end

%% Model the interest points as a series of bands on a linear probe object

[...
    ~, image_lengths, model_from_image_axes, model_from_image...
] = bilateralModel(...
    interest_points, point_alignment_outlier_threshold, true...
);

if display_model_from_image
    fg = figure; %#ok<UNRCH>
    imshow(I);
    plotBilateralModel( model_from_image, model_from_image_axes, [image_height, image_width], [], fg);
    title('Classified interest points (blue, black = tips; green, red = above/below first PCA component; yellow = unmatched)');
end

% The probe tips are marked with single points. All other probe segments
% must be marked with a point along each side of the probe.
%
% Since we are performing calibration, failure to satisfy this condition is
% an error. In a detection scenario, we would just eliminate the outlier
% points and continue.
if ~isempty(model_from_image.unmatched)
    error(sprintf(['Not all probe colour band junctions are marked with two points - One on each edge of the probe.\n',...
        'Consider increasing the outlier detection threshold used when pairing points.'])); %#ok<SPERR>
end

%% Match model extracted from the image to the user-supplied measurements of the probe
load(model_filename, 'probe');
if ~exist('probe', 'var')
    error('No variable called ''probe'' is loaded (which would contain probe measurements).')
end
    
image_to_measured_matches = matchProbeLengths(...
  probe.lengths, image_lengths, subject_gap_cost, query_gap_cost,...
  affine_weight, verbose_point_sequence_matching...
);

% Validate the matching, and flip the model extracted from the image if necessary
if any(~image_to_measured_matches)
    error('Some interest points in the image were not matched to known probe measurements.')
else
    image_to_measured_matches_expected = (1:length(image_to_measured_matches)).';
    aligned_forward = all(image_to_measured_matches == image_to_measured_matches_expected);
    aligned_reverse = all(image_to_measured_matches == flipud(image_to_measured_matches_expected));
    if ~xor(aligned_forward, aligned_reverse)
        error('Failed to find an ordered one-to-one mapping between interest points in the image and known probe measurements.')
    elseif aligned_reverse
        model_from_image_new = model_from_image;
        if isfield(model_from_image, 'head')
            model_from_image_new.tail = model_from_image.head;
        else
            model_from_image_new = rmfield(model_from_image_new, 'tail');
        end
        if isfield(model_from_image, 'tail')
            model_from_image_new.head = model_from_image.tail;
        else
            model_from_image_new = rmfield(model_from_image_new, 'head');
        end
        model_from_image_new.above = flipud(model_from_image_new.above);
        model_from_image_new.below = flipud(model_from_image_new.below);
        model_from_image = model_from_image_new;
    end
end

%% Locate regions corresponding to probe bands and probe colours

% Express the model in pixel coordinates as polygonal sections
n_bands = length(probe.lengths) - 1;
above = model_from_image.above;
below = model_from_image.below;
model_from_image_polygons = cell(n_bands, 1);
start_offset = 0;
if isfield(model_from_image, 'head')
    head = model_from_image.head;
    model_from_image_polygons{1} = [
            head(1:2);
            above(1, 1:2);
            below(1, 1:2)
        ];
    start_offset = 1;
end
for i = 1:(size(above, 1) - 1)
    model_from_image_polygons{i + start_offset} = [
            above(i:(i+1), 1:2);
            below((i+1):-1:i, 1:2)
        ];
end
if isfield(model_from_image, 'tail')
    tail = model_from_image.tail;
    model_from_image_polygons{end} = [
        above(end, 1:2);
        tail(1:2);
        below(end, 1:2);
    ];
end

% Obtain a mask for each band
probe_band_masks = false(image_height, image_width, n_bands);
for i = 1:n_bands
    probe_band_masks(:, :, i) = roipoly(...
            I, model_from_image_polygons{i}(:, 1), model_from_image_polygons{i}(:, 2)...
        );
end

if display_probe_band_masks
    for i = 1:n_bands %#ok<UNRCH>
        figure
        imshow(probe_band_masks(:, :, i));
        title(sprintf('Mask for probe band %d', i))
    end
end

% Group bands by colour
colors_filter = (probe.colors ~= 0);
probe_colors = unique(probe.colors(colors_filter));
if any(diff(probe_colors) ~= 1)
    error('Probe bands must be given colour labels that are consecutive integers.')
end
n_colors = length(probe_colors);
probe_color_masks = false(image_height, image_width, n_colors);
for i = 1:n_bands
    if colors_filter(i)
        probe_color_masks(:, :, probe.colors(i)) =...
            probe_color_masks(:, :, probe.colors(i)) | probe_band_masks(:, :, i);
    end
end

if display_probe_color_masks
    for i = 1:n_colors %#ok<UNRCH>
        figure
        imshow(probe_color_masks(:, :, i));
        title(sprintf('Mask for probe colour %d', i))
    end
end

%% Create photometric invariant representations of the probe colours

% Obtain hue values
H = rgb2hue(I);
    
if display_hue_image
    figure %#ok<UNRCH>
    H_color = ones(image_height, image_width, image_n_channels);
    H_color(:, :, 1) = H;
    H_color = hsv2rgb(H_color);
    imshowpair(H, H_color, 'montage');
    title('Hue channel of original image')
end

% Compute hue Gaussian density estimators from probe colors

probe_color_distributions = zeros(...
        probe_color_distribution_resolution, n_colors...
    );
for i = 1:n_colors
    [...
        probe_color_distributions(:, i),...
        probe_color_distribution_increment...
    ] = hueGaussianDensityEstimator(...
        H, probe_color_masks(:, :, i), probe_color_distribution_resolution...
    );
end

if plot_hue_estimators
    legend_names = cell(n_colors, 1); %#ok<UNRCH>
    for i = 1:n_colors
        legend_names{i} = sprintf('Probe color %d', i);
    end
    plotHueDensityEstimator(...
        probe_color_distribution_increment, probe_color_distributions, legend_names...
    );
    title('Hue Gaussian density estimators for colors on the probe')
end

%% Create a hue classifier

probe_color_classifier = mlDiscreteClassifier(...
    probe_color_distributions, probe_color_distribution_increment,...
    'periodic'...
);

if plot_hue_classifier
    plotHueClassifier(...
        probe_color_distribution_increment, probe_color_classifier,...
        n_colors...
    );
    title('Hue classifier for colors on the probe')
end

%% Save results to a file
save_variables_list = [ parameters_list, {...
        'probe',...
        'probe_color_distributions',...
        'probe_color_distribution_increment',...
        'probe_color_classifier'
    } ];
uisave(save_variables_list,'probeDetectionModel')
disp('Reminder: The output model is specific to the probe, camera, and camera parameters.')