%% RGB Noise Estimation
% Create a rough model of RGB standard deviations from one or more images.
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
% ### Probe detection model
% A '.mat' file containing the following variable:
% - 'rgb_sigma_polyfit': An n x n_channels array, containing the coefficients
%   of polynomials fitted to `(mean(channel_values),std(channel_values))`
%   pairs. 'n' is the degree of the polynomials, and 'n_channels' is 3 (for
%   Red, Green, Blue, in that order). 'rgb_sigma_polyfit' describes the
%   relationship between color channel values and color channel noise, as
%   calculated using MATLAB's `polyfit` function. Specifically, each column
%   of 'rgb_sigma_polyfit' corresponds to the first output argument, 'p',
%   of `polyfit`.
%
% ## References
% - T. Gevers and H. Stokman. "Robust Histogram Construction from Color
%   Invariants for Object Recognition". IEEE Transactions on Pattern
%   Analysis and Machine Intelligence, vol. 26, no. 1, pp. 113-118, Jan.
%   2004.

% Bernard Llanos
% Spring 2016 research assistantship supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created August 9, 2016

%% Initialization
if ~exist('choice', 'var') || ~strcmpi(choice, 'z')
    choice = 'y';
    mu = zeros(0, 3);
    sigma = zeros(0, 3);
    fg = [];
    n_channels = 3;
    channels = cell(n_channels, 1);
    polyfit_degree_default = 1;
    polyfit_degree = input(sprintf(...
            'Desired degree of fitted polynomial (Default %d): ',...
            polyfit_degree_default)...
        );
    if isempty(polyfit_degree)
        polyfit_degree = polyfit_degree_default;
    end
    n_points_desired = polyfit_degree + 1;
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
        I = imread(fullfile(pathname, filename));
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
            mask = roipoly;
            if ~isempty(mask)
                mu(end + 1, :) = zeros(1, 3); %#ok<SAGROW>
                sigma(end + 1, :) = zeros(1, 3); %#ok<SAGROW>
                for i = 1:n_channels
                    px = channels{i}(mask);
                    mu(end, i) = mean(px);
                    sigma(end, i) = std(px);
                end
            else
                % 'Resume after adjusting view' allows the user to zoom into the
                % image, then re-run the script to continue data collection.
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
    rgb_sigma_polyfit = zeros(polyfit_degree + 1, n_channels);
    for i = 1:n_channels
        rgb_sigma_polyfit(:, i) = polyfit(mu(:, i), sigma(:, i), polyfit_degree);
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
        scatter(mu(:, i), sigma(:, i), markerspec{i});
        plot(x, polyval(rgb_sigma_polyfit(:, i), x), linespec{i});
    end
    hold off
    title('Region RGB standard deviations')
    xlabel('Mean colour value')
    ylabel('Standard deviation')
    
    % Save to a file
    uisave('rgb_sigma_polyfit','rgbStddev')
end