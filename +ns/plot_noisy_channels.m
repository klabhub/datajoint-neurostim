function plot_noisy_channels(noisyChannels, signal, timepoints, varargin)
% PLOT_NOISYCHANNELS Visualizes results from find_noisy_channels.
%
%   plot_noisyChannels(noisyChannels, signal, timepoints, Name, Value, ...)
%
%   Creates two figures:
%   1. Bad Channels by Category: Shows plots of channels flagged by each
%      criterion in separate vertical tiles. Channels are plotted overlaid.
%   2. All vs Good Channels: Shows a plot of all channels (overlaid) in the
%      top tile and only 'good' channels (overlaid) in the bottom tile.
%
%   Inputs:
%       noisyChannels - Structure output from find_noisy_channels function.
%                       Must contain fields like .all, .badBy...
%       signal        - (nChannels x nTimepoints) EEG data matrix used in
%                       find_noisy_channels.
%       timepoints    - (1 x nTimepoints) Time vector corresponding to signal.
%       varargin      - Optional Name-Value pairs:
%           'TimeRange'        - [startTime, endTime] vector to plot only a
%                                specific time segment. Default: plots all.
%           'MaxChannelsPerPlot' - Maximum number of channels to display in
%                                  a single tile (subset chosen if exceeded).
%                                  Default: 50.
%           'YLimIQRScale'     - Factor to multiply IQR by for automatic Y limits
%                                (Median +/- Factor * IQR). Default: 2.0.
%                                Ignored if 'YLim' is provided.
%           'YLim'             - [minY, maxY] vector to manually set Y-axis
%                                limits for ALL tiles. Overrides automatic scaling.
%                                Default: [].
%           'LineWidth'        - Line width for plots. Default: 1.
%           'LineAlpha'        - Transparency level for plotted lines (0=clear, 1=opaque).
%                                Default: 1.0.
%           'ChannelLocations' - (nChannels x 3)
%           'FigureVisible'    - 'on' (default) or 'off' to show figures.
%
%   Requires: MATLAB R2019b or later for tiledlayout.

% --- Input Parser ---
p = inputParser;
p.FunctionName = 'plot_noisy_channels';

% Validation functions
isValidNoisyStruct = @(x) isstruct(x) && isfield(x, 'all') && ~isempty(regexp(strjoin(fieldnames(x)), 'badBy', 'once'));
isValidSignal = @(x) isnumeric(x) && ismatrix(x) && ~isempty(x);
isValidTimepoints = @(x, sig) isnumeric(x) && isvector(x) && length(x) == size(sig, 2);
isValidRange = @(x) isempty(x) || (isnumeric(x) && numel(x) == 2 && x(1) < x(2));
isValidMaxChan = @(x) isnumeric(x) && isscalar(x) && x > 0 && floor(x)==x;
isValidIQRScale = @(x) isnumeric(x) && isscalar(x) && x > 0;
isValidYLim = @(x) isempty(x) || (isnumeric(x) && isvector(x) && numel(x) == 2 && x(1) < x(2)); % Validation for YLim
isValidLineWidth = @(x) isnumeric(x) && isscalar(x) && x > 0;
isValidAlpha = @(x) isnumeric(x) && isscalar(x) && x >= 0 && x <= 1; % Validation for LineAlpha
isValidVisibility = @(x) ischar(x) || isstring(x) && ismember(lower(x), {'on', 'off'});
isValidChannelLocations = @(x) isempty(x) || (ismatrix(x) && size(x,2)==3);

% Add inputs
addRequired(p, 'noisyChannels', isValidNoisyStruct);
addRequired(p, 'signal', isValidSignal);
addRequired(p, 'timepoints', @(x) isValidTimepoints(x, signal)); % Pass signal for validation
addParameter(p, 'TimeRange', [], isValidRange);
addParameter(p, 'MaxChannelsPerPlot', 50, isValidMaxChan);
addParameter(p, 'YLimIQRScale', 2.0, isValidIQRScale);
addParameter(p, 'YLim', [], isValidYLim); % New YLim parameter
addParameter(p, 'LineWidth', 1, isValidLineWidth);
addParameter(p, 'LineAlpha', 1.0, isValidAlpha); % New LineAlpha parameter
addParameter(p, 'FigureVisible', 'on', isValidVisibility);
addParameter(p, 'ChannelLocations', [], isValidChannelLocations);

% Parse
parse(p, noisyChannels, signal, timepoints, varargin{:});
params = p.Results;

% --- Data Preparation ---
[nChannels, nTimepoints] = size(signal);
% Fs = 1 / median(diff(timepoints)); % Estimate sampling frequency (not strictly needed here)

% Time selection
if isempty(params.TimeRange)
    timeIdx = 1:nTimepoints;
    timeToPlot = timepoints;
    timeTitleSuffix = '';
else
    timeIdx = find(timepoints >= params.TimeRange(1) & timepoints <= params.TimeRange(2));
    if isempty(timeIdx)
        warning('plot_noisyChannels:TimeRange', 'Specified TimeRange [%.2f, %.2f] is outside data bounds or yields no points. Plotting all time.', params.TimeRange(1), params.TimeRange(2));
        timeIdx = 1:nTimepoints;
        timeToPlot = timepoints;
        timeTitleSuffix = '';
    else
        timeToPlot = timepoints(timeIdx);
        timeTitleSuffix = sprintf(' (Time: %.2f-%.2f s)', params.TimeRange(1), params.TimeRange(2));
    end
end

% --- Figure 1: Bad Channels by Category ---
fprintf('Generating Figure 1: Bad Channels by Category (Overlaid)...\n');

fig1 = figure('Name', 'Bad Channels by Category (Overlaid)', 'Visible', params.FigureVisible, 'NumberTitle', 'off');
drawnow; % Ensure figure exists before adding tiles

% Identify categories and channels
fields = fieldnames(noisyChannels);
badCategories = {};
badIndicesByCategory = {};
for i = 1:length(fields)
    if startsWith(fields{i}, 'badBy') && ~strcmp(fields{i}, 'all')
        indices = noisyChannels.(fields{i});
        if ~isempty(indices)
            indices = indices(indices >= 1 & indices <= nChannels);
             if ~isempty(indices)
                badCategories{end+1} = fields{i}; %#ok<AGROW>
                badIndicesByCategory{end+1} = unique(indices(:)'); %#ok<AGROW>
             end
        end
    end
end

nCategories = length(badCategories);

if nCategories == 0
    clf(fig1);
    axes(fig1);
    text(0.5, 0.5, 'No channels found in specific bad categories.', 'HorizontalAlignment', 'center', 'FontSize', 12);
    title('Bad Channels by Category');
    set(gca, 'XTick', [], 'YTick', []);
else
    t1 = tiledlayout(fig1, nCategories, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

    for iCat = 1:nCategories
        categoryName = badCategories{iCat};
        categoryIndices = badIndicesByCategory{iCat};
        nCatChans = length(categoryIndices);

        % Select subset if too many channels
        plotIndices = categoryIndices;
        titleSuffix = '';
        if nCatChans > params.MaxChannelsPerPlot
            plotIndices = categoryIndices(1:params.MaxChannelsPerPlot);
            titleSuffix = sprintf(' (showing %d of %d)', params.MaxChannelsPerPlot, nCatChans);
            warning('plot_noisyChannels:MaxChan', 'Category "%s": Too many channels (%d). Plotting only the first %d.', categoryName, nCatChans, params.MaxChannelsPerPlot);
        end

        nToPlot = length(plotIndices);
        signalToPlot = signal(plotIndices, timeIdx);

        % Plotting (No Offset)
        ax = nexttile(t1);
        hold(ax, 'on');
        % Plot channels directly on top of each other
        hLines = plot(ax, timeToPlot, signalToPlot', 'LineWidth', params.LineWidth);
        
        % Apply Transparency
        if params.LineAlpha < 1.0 && ~isempty(hLines)
            applyTransparency(ax, hLines, params.LineAlpha);
        end
        
        hold(ax, 'off');

        % Set Axes and Labels
        axis(ax, 'tight'); % Tighten X axis first
        xlabel(ax, 'Time (s)');
        ylabel(ax, 'Amplitude (\muV)'); % Standard Y label
        title(ax, sprintf('%s (%d Channels)%s', formatCategoryName(categoryName), nCatChans, titleSuffix));
        grid(ax, 'on');

        % Set Y Limits (Manual or Robust)
        if ~isempty(params.YLim)
            ylim(ax, params.YLim); % Apply manual limits
        else
            % Only apply robust limits if manual ones aren't given
            setRobustYLim(ax, signalToPlot, params.YLimIQRScale);
        end
    end
    title(t1, ['Bad Channels by Category' timeTitleSuffix], 'FontWeight', 'bold');
end
drawnow;

% --- Figure 2: All vs Good Channels (Overlaid) ---
fprintf('Generating Figure 2: All vs Good Channels (Overlaid)...\n');

fig2 = figure('Name', 'All vs Good Channels (Overlaid)', 'Visible', params.FigureVisible, 'NumberTitle', 'off');
drawnow;

t2 = tiledlayout(fig2, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

% Tile 1: All Channels
ax1 = nexttile(t2);
fprintf(' Plotting all channels (overlaid)...\n');
allIndices = 1:nChannels;
nAllChans = nChannels;
plotIndicesAll = allIndices;
titleSuffixAll = '';
if nAllChans > params.MaxChannelsPerPlot
    plotIndicesAll = sort(randperm(nAllChans, params.MaxChannelsPerPlot)); % Plot random subset
    titleSuffixAll = sprintf(' (showing random %d of %d)', params.MaxChannelsPerPlot, nAllChans);
     warning('plot_noisyChannels:MaxChan', 'All Channels: Too many channels (%d). Plotting only a random subset of %d.', nAllChans, params.MaxChannelsPerPlot);
end
nToPlotAll = length(plotIndicesAll);
signalToPlotAll = signal(plotIndicesAll, timeIdx);

% Plotting (No Offset)
hold(ax1, 'on');
hLines1 = plot(ax1, timeToPlot, signalToPlotAll', 'LineWidth', params.LineWidth);
% Apply Transparency
if params.LineAlpha < 1.0 && ~isempty(hLines1)
    applyTransparency(ax1, hLines1, params.LineAlpha);
end
hold(ax1, 'off');

% Set Axes and Labels
axis(ax1, 'tight');
xlabel(ax1, 'Time (s)');
ylabel(ax1, 'Amplitude (\muV)');
title(ax1, sprintf('All Channels (%d)%s', nAllChans, titleSuffixAll));
grid(ax1, 'on');

% Set Y Limits (Manual or Robust)
if ~isempty(params.YLim)
    ylim(ax1, params.YLim); % Apply manual limits
else
    setRobustYLim(ax1, signalToPlotAll, params.YLimIQRScale);
end


% Tile 2: Good Channels Only
ax2 = nexttile(t2);
fprintf(' Plotting good channels (overlaid)...\n');
goodIndices = setdiff(1:nChannels, noisyChannels.all);
nGoodChans = length(goodIndices);

if nGoodChans == 0
    text(ax2, 0.5, 0.5, 'No good channels found.', 'HorizontalAlignment', 'center', 'FontSize', 12);
    title(ax2, 'Good Channels Only (0)');
    set(ax2, 'XTick', [], 'YTick', []);
else
    plotIndicesGood = goodIndices;
    titleSuffixGood = '';
    if nGoodChans > params.MaxChannelsPerPlot
        plotIndicesGood = goodIndices(sort(randperm(nGoodChans, params.MaxChannelsPerPlot))); % Plot random subset
        titleSuffixGood = sprintf(' (showing random %d of %d)', params.MaxChannelsPerPlot, nGoodChans);
         warning('plot_noisyChannels:MaxChan', 'Good Channels: Too many channels (%d). Plotting only a random subset of %d.', nGoodChans, params.MaxChannelsPerPlot);
    end
    nToPlotGood = length(plotIndicesGood);
    signalToPlotGood = signal(plotIndicesGood, timeIdx);

    % Plotting (No Offset)
    hold(ax2, 'on');
    hLines2 = plot(ax2, timeToPlot, signalToPlotGood', 'LineWidth', params.LineWidth);
    % Apply Transparency
    if params.LineAlpha < 1.0 && ~isempty(hLines2)
        applyTransparency(ax2, hLines2, params.LineAlpha);
    end
    hold(ax2, 'off');

    % Set Axes and Labels
    axis(ax2, 'tight');
    xlabel(ax2, 'Time (s)');
    ylabel(ax2, 'Amplitude (\muV)');
    title(ax2, sprintf('Good Channels Only (%d)%s', nGoodChans, titleSuffixGood));
    grid(ax2, 'on');

    % Set Y Limits (Manual or Robust)
    if ~isempty(params.YLim)
        ylim(ax2, params.YLim); % Apply manual limits
    else
        setRobustYLim(ax2, signalToPlotGood, params.YLimIQRScale);
    end

    % Link time axes
    linkaxes([ax1, ax2], 'x');

    % topoplot
    if ~isempty(params.ChannelLocations)
        % scale to unit sphere
        % rad = sqrt(sum(params.ChannelLocations.^2,2));
        % params.ChannelLocations = params.ChannelLocations ./ rad;
        plot_topography_with_noise(params.ChannelLocations, noisyChannels);
    end
    
end

title(t2, ['All Channels vs Good Channels Comparison (Overlaid)' timeTitleSuffix], 'FontWeight', 'bold');
drawnow;

fprintf('Plotting finished.\n');

end % End of main function

% --- Helper Function: Format Category Name ---
function formattedName = formatCategoryName(rawName)
    % Converts 'badByHFNoise' to 'Bad by HF Noise' etc.
    formattedName = regexprep(rawName, 'badBy(.)(.*)', 'Bad by ${upper($1)}$2');
    formattedName = regexprep(formattedName, '([A-Z])', ' $1');
    formattedName = strrep(formattedName, ' H F ', ' HF ');
    formattedName = strrep(formattedName, ' S N R', ' SNR');
    formattedName = strrep(formattedName, ' Ransac', ' RANSAC');
    formattedName = strtrim(formattedName);
end

% --- Helper Function: Set Robust Y Limits ---
function setRobustYLim(ax, dataToPlot, iqrScale)
    % Sets the Y limits based on Median +/- iqrScale * IQR of the data.
    
    allDataInTile = dataToPlot(:);
    allDataInTile = allDataInTile(isfinite(allDataInTile));

    if ~isempty(allDataInTile)
        medVal = median(allDataInTile);
        iqrVal = iqr(allDataInTile);

        if iqrVal < 1e-6 % Handle flat data
             yLower = medVal - 1; 
             yUpper = medVal + 1;
        else
            yLower = medVal - iqrScale * iqrVal;
            yUpper = medVal + iqrScale * iqrVal;
        end

         % Add padding if range is still very small
         if abs(yUpper - yLower) < 1e-5
             padding = max(0.5, abs(medVal)*0.05 + 0.1); % Add padding
             yLower = yLower - padding;
             yUpper = yUpper + padding;
         end
         % Ensure lower < upper
         if yLower >= yUpper
              yLower = yUpper - 0.1; % Ensure separation
         end

        ylim(ax, [yLower, yUpper]);
    % else leave ylim as determined by 'axis tight'
    end
end

% --- Helper Function: Apply Transparency ---
function applyTransparency(ax, hLines, alphaValue)
    % Applies transparency to lines plotted on axes.
    % Uses the axes ColorOrder to determine original colors.
    if alphaValue >= 1.0 || alphaValue < 0 || isempty(hLines)
        return; % No transparency needed or possible
    end
    
    currentColors = get(ax, 'ColorOrder');
    numColors = size(currentColors, 1);
    
    for iLine = 1:length(hLines)
         try
             % Get the color index used for this line
             colorIdx = mod(iLine-1, numColors) + 1; 
             lineColor = currentColors(colorIdx, :); 
             % Set the Color property using RGBA format
             hLines(iLine).Color = [lineColor, alphaValue];
         catch ME
              warning('plot_noisyChannels:TransparencyError', ...
                      'Could not apply transparency to line %d. Error: %s', iLine, ME.message);
              % Don't apply transparency if it fails for some reason
         end
    end
end

function figHandle = plot_topography_with_noise(chanLocs, noisyChannels, varargin)
% PLOT_TOPOGRAPHY_WITH_NOISE Displays channel locations with noise categories.

figHandle = gcf;

ch_list = noisyChannels;
ch_cats = string(fieldnames(ch_list));
for catN = ch_cats'

    if ~startsWith(catN, "badBy") || isempty(ch_list.(catN))

        ch_list.(catN) = [];

    end

end
ephys.egi.topoplot3d(chanLocs, "annotate", ch_list);

fprintf('Topoplot created.\n');

end % END OF FUNCTION plot_topography_with_noise