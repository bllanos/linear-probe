%% RGB Noise Estimation (Interactive)
% Create a rough model of RGB standard deviations from one or more images,
% by interactively selecting regions of constant colour.
%
% ## Usage
%   Run the script and follow the prompts.
%
% ## Input
% The script will prompt for the locations of one (or, optionally, several)
% images to use for estimating RGB noise. The images should have been taken
% under the same camera parameters and lighting conditions. For accurate
% noise estimation, mark regions in the images that correspond to
% constant-coloured surfaces with little variation in shading/highlighting.
%
% ## Output
%
% ### RGB noise model
% A '.mat' file containing the following variables:
% - 'rgb_sigma_polyfit': A 2 x n_channels array, containing the coefficients
%   of lines fitted to `(mean(channel_values),var(channel_values))` pairs.
%   'n_channels' is 3 (for Red, Green, Blue, in that order).
%   'rgb_sigma_polyfit' describes the relationship between colour channel
%   values and colour channel noise, as calculated using MATLAB's `polyfit`
%   function. Specifically, each column of 'rgb_sigma_polyfit' corresponds
%   to the first output argument, 'p', of `polyfit`.
%
%   Note that colour channel values are in the range [0, 1] (and standard
%   deviations correspond to this range as well, therefore).
% - 'image_filenames': A cell vector containing the full filepaths of the
%   images used to compute the RGB noise model, for reference.
% - 'image_roi_polygons': A cell vector with the same length as
%   'image_filenames', where the i-th element is a cell vector of k x 2 arrays
%   representing the polygonal regions of interest used to calculate RGB
%   noise datapoints from the image `image_filenames{i}`. 'k' is the number
%   of vertices in the polygon, and the two columns store pixel x and y
%   coordinates, respectively.
%
% ## References
% - T. Gevers and H. Stokman. "Robust Histogram Construction from Color
%   Invariants for Object Recognition". IEEE Transactions on Pattern
%   Analysis and Machine Intelligence, vol. 26, no. 1, pp. 113-118, Jan.
%   2004.
% - Martinec, E. (2008). Noise, dynamic range and bit depth in digital SLR.
%   Retrieved from http://theory.uchicago.edu/∼ejm/pix/20d/tests/noise/

% Bernard Llanos
% Spring 2016 research assistantship supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created August 9, 2016

%% Initialization
if ~exist('choice', 'var') || ~strcmpi(choice, 'z')
    choice = 'y';
    mu = zeros(0, 3);
    px_var = zeros(0, 3);
    fg = [];
    n_channels = 3;
    channels = cell(n_channels, 1);
    n_points_desired = 2;
    image_filenames = {};
    image_roi_polygons = {};
else
    disp('Enter the next region of interest on the figure.')
    choice = [];
end

%% Collect data from image regions of interest
while size(mu, 1) < n_points_desired || ~strcmpi(choice, 'n')
    if strcmpi(choice, 'y')
        if ishandle(fg)
            close(fg)
            fg = [];
        end
        [filename, pathname] = uigetfile(...
            {'*.png;*.jpg;*.bmp;*.tif;*.gif','All Image Files'; '*.*','All Files' },...
            'Select an image file'...
            );
        if ~filename
            break
        end
        image_filenames{end + 1} = fullfile(pathname, filename); %#ok<SAGROW>
        I = imread(image_filenames{end});
        if size(I, 3) ~= n_channels
            error('Expected a 3-channel RGB image.')
        else
            J = im2double(I);
            for i = 1:n_channels
                channels{i} = J(:, :, i);
            end
            fg = figure;
            imshow(I);
            choice = [];
        end
        
    elseif isempty(choice)
        try
            figure(fg);
            title(sprintf('Select region %d (of at least %d) of constant colour (Esc, etc., to cancel)',...
            size(mu, 1) + 1, n_points_desired))
            [mask, x_poly, y_poly] = roipoly;
            if ~isempty(mask)
                mu(end + 1, :) = zeros(1, 3); %#ok<SAGROW>
                px_var(end + 1, :) = zeros(1, 3); %#ok<SAGROW>
                for i = 1:n_channels
                    px = channels{i}(mask);
                    mu(end, i) = mean(px);
                    px_var(end, i) = var(px);
                end
                if length(image_roi_polygons) < length(image_filenames)
                    image_roi_polygons{length(image_filenames)} = {[x_poly, y_poly]}; %#ok<SAGROW>
                else
                    image_roi_polygons{end}{end + 1} = [x_poly, y_poly];
                end
            else
                % 'Resume after adjusting view' allows the user to zoom into the
                % image, then re-run the script to continue data collection.
                %
                % While more complicated, a GUI to display the image and
                % offer buttons for adjusting the region of the image being
                % displayed and for marking regions of interest would be
                % ideal for this purpose.
                choice = input('New image? (y), Resume after adjusting view (z), or Stop data collection? (n): ', 's');
                if strcmpi(choice, 'z')
                    break
                end
            end
        
        catch ME
            if strcmp(ME.identifier, 'MATLAB:hgbuiltins:object_creation:InvalidConvenienceArgHandle')
                choice = []; %#ok<NASGU>
                error('Figure is no longer valid. Cannot resume data collection.');
            end
        end
        
    else
        fprintf('Unrecognized choice "%s"\n', choice);
        break;
    end
end

%% Process results
if strcmpi(choice, 'z')
    disp('Re-run the script after adjusting the figure.')
elseif size(mu, 1) < n_points_desired
    disp('Insufficient data collected.')
    if ishandle(fg)
        close(fg);
    end
else
    
    % Fit curves to the data
    rgb_sigma_polyfit = zeros(2, n_channels);
    for i = 1:n_channels
        rgb_sigma_polyfit(:, i) = polyfit(mu(:, i), px_var(:, i), 1);
    end
    
    % Plot the results
    if ~ishandle(fg)
        fg = figure;
    else
        clf(fg)
    end
    hold on
    markerspec = {'ro', 'go', 'bo'};
    linespec = {'r--', 'g--', 'b--'};
    x = linspace(0, 1, 50);
    for i = 1:n_channels
        scatter(mu(:, i), sqrt(px_var(:, i)), markerspec{i});
        plot(x, sqrt(polyval(rgb_sigma_polyfit(:, i), x)), linespec{i});
    end
    hold off
    title('Region RGB standard deviations')
    xlabel('Mean colour value')
    ylabel('Standard deviation')
    
    % Save to a file
    uisave({'rgb_sigma_polyfit', 'image_filenames', 'image_roi_polygons'},'rgbStddev')
end