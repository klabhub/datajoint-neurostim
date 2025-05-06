function H = topoplot3d(chanLocs, options)
%TOPLOT3D Plots 3D channel locations with optional surface and annotation.
%
%   H = TOPLOT3D(chanLocs) plots the 3D electrode locations specified in
%   the n x 3 matrix chanLocs (X, Y, Z coordinates). It returns a struct H
%   containing handles to the graphical objects. A default surface based
%   on the convex hull is drawn.
%
%   H = TOPLOT3D(chanLocs, Name, Value, ...) allows specifying optional
%   parameters via Name-Value pairs.
%
%   Optional Name-Value Pairs:
%
%   'SurfaceStyle':      Style of the surface ('trisurf' or 'off').
%                        'trisurf' (default) draws a surface based on the
%                        convex hull of the points. 'off' draws no surface.
%   'SurfaceColor':      RGB triplet or color character (e.g., 'g', [0.5 0.5 0.5])
%                        for the surface face color. Default is light gray.
%   'SurfaceAlpha':      Transparency of the surface (0 to 1). Default is 0.6.
%   'MarkerSize':        Size of the markers for electrode locations. Default is 30.
%   'DefaultMarkerColor': Color for markers not part of an annotation group.
%                         Default is 'k' (black).
%   'annotate':  A struct used for grouping and coloring channels.
%                        - Field names are arbitrary group names (e.g., 'Frontal', 'Parietal').
%                        - Field values are cell arrays of channel label strings
%                          belonging to that group (e.g., {'Fp1', 'Fz', 'Fp2'}).
%                        Channels within the same group are plotted with the
%                        same marker color and annotated with text labels.
%                        Default is an empty struct (no annotation).
%   'ChannelLabels':     A cell array of strings containing the label for *each*
%                        channel in chanLocs, in the *same order*. This is
%                        **REQUIRED** if 'AnnotationStruct' is used. Default is {}.
%   'AnnotationFontSize': Font size for the text labels used in annotation. Default is 8.
%   'LabelOffset':       Small offset applied to text labels to avoid overlap
%                        with markers. Default is 0.02.
%   'AxesHandle':        Handle to an existing axes object where the plot
%                        should be drawn. If empty, a new figure and axes
%                        are created. Default is [].
%
%   Output:
%   H: A struct containing handles to the plotted elements:
%      H.Surface:       Handle to the surface object (if plotted).
%      H.AllPoints:     Handle to the scatter plot of all points (if not annotating).
%      H.DefaultPoints: Handle to scatter plot of points *not* in any annotation group.
%      H.GroupPoints:   Struct with fields corresponding to annotation group names.
%                       Each field contains the handle to the scatter plot for that group.
%      H.GroupLabels:   Struct with fields corresponding to annotation group names.
%                       Each field contains an array of handles to the text labels for that group.
%      H.Axes:          Handle to the axes where the plot was drawn.
%
%   Example Usage:
%       % 1. Basic plot
%       locs = rand(64, 3) * 10; % Example locations
%       figure;
%       H1 = topoplot3d(locs);
%
%       % 2. Plot without surface, different marker size
%       figure;
%       H2 = topoplot3d(locs, 'SurfaceStyle', 'off', 'MarkerSize', 50);
%
%       % 3. Plot with annotation
%       labels = cell(64, 1); % Create matching labels
%       for i = 1:64, labels{i} = sprintf('Ch%d', i); end
%
%       annotationGroups = struct();
%       annotationGroups.Frontal = {'Ch1', 'Ch2', 'Ch3', 'Ch4'};
%       annotationGroups.Occipital = {'Ch60', 'Ch61', 'Ch62', 'Ch63', 'Ch64'};
%
%       figure;
%       H3 = topoplot3d(locs, 'AnnotationStruct', annotationGroups, ...
%                           'ChannelLabels', labels, 'SurfaceAlpha', 0.3);
%
%   Requires MATLAB R2019b or later for the 'arguments' block.

% --- Input Argument Parsing and Validation ---
arguments
    chanLocs (:,3) {mustBeNumeric, mustBeReal}

    % --- Options ---
    options.SurfaceStyle (1,1) string {mustBeMember(options.SurfaceStyle, ["trisurf", "off"])} = "trisurf"
    options.SurfaceColor {validatecolor} = .99*[1 1 1] % Default light gray
    options.SurfaceAlpha (1,1) {mustBeNumeric, mustBeReal, mustBeInRange(options.SurfaceAlpha, 0, 1)} = 0.9
    options.MarkerSize (1,1) {mustBeNumeric, mustBePositive} = 50
    options.DefaultMarkerColor {validatecolor} = 'k' % Default white
    options.annotate (1,1) struct = struct()
    options.ChannelLabels (:,1) cell = {}
    options.AnnotationFontSize (1,1) {mustBeNumeric, mustBePositive} = 8
    options.LabelOffset (1,1) {mustBeNumeric, mustBeReal} = 0.02
    options.AxesHandle = gobjects(0) % Use gobjects(0) as default for optional handle
end

% --- Initial Setup ---
H = struct(); % Initialize handles structure

% Determine target axes
if isempty(options.AxesHandle) || ~isgraphics(options.AxesHandle, 'axes')
    figure;
    ax = gca;
else
    ax = options.AxesHandle;
end
H.Axes = ax; % Store axes handle

% Extract coordinates for convenience
x = chanLocs(:, 1);
y = chanLocs(:, 2);
z = chanLocs(:, 3);
nChannels = size(chanLocs, 1);

% --- Plot Surface (Optional) ---
hold(ax, 'on'); % Hold axes for subsequent plots

if options.SurfaceStyle == "trisurf"
    try
        K = convhull(x, y, z); % Calculate convex hull triangulation
        h_surf = trisurf(K, x, y, z, 'Parent', ax); % Plot triangulated surface

        % Apply surface styling
        set(h_surf, 'FaceColor', options.SurfaceColor, ...
                    'EdgeColor', 'none', ...
                    'FaceAlpha', options.SurfaceAlpha);
        H.Surface = h_surf; % Store handle
    catch ME
        warning('Could not compute or plot convex hull surface: %s', ME.message);
        H.Surface = gobjects(0);
    end
else
    H.Surface = gobjects(0); % Indicate no surface handle
end

% --- Annotation Setup ---
doAnnotation = ~isempty(fieldnames(options.annotate));
isAnnotated = false(nChannels, 1); % Keep track of channels handled by annotation

if doAnnotation
    % Validation: ChannelLabels must be provided and match chanLocs
    if isempty(options.ChannelLabels)
        options.ChannelLabels = string(1:nChannels);
    elseif numel(options.ChannelLabels) ~= nChannels
        error('Number of ChannelLabels must match the number of rows in chanLocs.');
    end

    % Create a map from label string to index for quick lookup
    labelMap = containers.Map(options.ChannelLabels, 1:nChannels);

    % Get group names and assign colors
    groupNames = fieldnames(options.annotate);
    nGroups = numel(groupNames);
    % Use 'lines' colormap for distinct colors, repeat if more groups than colors
    groupColors = lines(nGroups);

    % Initialize handles substructs
    H.GroupPoints = struct();
    H.GroupLabels = struct();
end

% --- Plot Points and Annotate ---

if doAnnotation
    fprintf('Annotating channel groups...\n');
    for iGroup = 1:nGroups
        groupName = groupNames{iGroup};
        groupColor = groupColors(iGroup, :);
        groupLabels = string(options.annotate.(groupName));

        if ~iscell(groupLabels) % Ensure it's a cell array
             groupLabels = cellstr(groupLabels);
        end

        % Find indices corresponding to labels in this group
        idxGroup = [];
        validLabelsInGroup = {};
        for iLabel = 1:numel(groupLabels)
            currentLabel = groupLabels{iLabel};
            if isKey(labelMap, currentLabel)
                idx = labelMap(currentLabel);
                % Avoid adding index if already annotated by a previous group
                if ~isAnnotated(idx)
                    idxGroup = [idxGroup; idx]; %#ok<AGROW>
                    validLabelsInGroup = [validLabelsInGroup; currentLabel]; %#ok<AGROW>
                    isAnnotated(idx) = true; % Mark as handled
                else
                     fprintf('  Warning: Label "%s" in group "%s" was already assigned to a previous group. Skipping.\n', currentLabel, groupName);
                end
            else
                fprintf('  Warning: Label "%s" in group "%s" not found in ChannelLabels. Skipping.\n', currentLabel, groupName);
            end
        end

        if isempty(idxGroup)
            fprintf('  No valid, unassigned channels found for group "%s".\n', groupName);
            H.GroupPoints.(groupName) = gobjects(0);
            H.GroupLabels.(groupName) = gobjects(0);
            continue; % Skip to next group if no valid channels found
        end

        fprintf('  Plotting group "%s" (%d channels) with color [%.2f %.2f %.2f]\n', ...
                groupName, numel(idxGroup), groupColor(1), groupColor(2), groupColor(3));

        % Plot points for this group
        h_pts = scatter3(ax, x(idxGroup), y(idxGroup), z(idxGroup), ...
                         options.MarkerSize, groupColor, 'filled');
        H.GroupPoints.(groupName) = h_pts;

        % Add text labels for this group
        h_labels_group = gobjects(numel(idxGroup), 1);
        for i = 1:numel(idxGroup)
            idx = idxGroup(i);
            labelStr = validLabelsInGroup{i}; % Use the validated label
            h_labels_group(i) = text(ax, x(idx) + options.LabelOffset, ...
                                     y(idx) + options.LabelOffset, ...
                                     z(idx) + options.LabelOffset, ...
                                     labelStr, ...
                                     'FontSize', options.AnnotationFontSize, ...
                                     'Color', groupColor, ...
                                     'Interpreter', 'none'); % Avoid TeX interpretation
        end
        H.GroupLabels.(groupName) = h_labels_group;
    end

    % Plot remaining points (not annotated) with default color
    idxDefault = find(~isAnnotated);
    if ~isempty(idxDefault)
        fprintf('Plotting %d remaining channels with default color.\n', numel(idxDefault));
        h_def_pts = scatter3(ax, x(idxDefault), y(idxDefault), z(idxDefault), ...
                             options.MarkerSize, options.DefaultMarkerColor, 'filled');
        H.DefaultPoints = h_def_pts;
    else
         H.DefaultPoints = gobjects(0);
    end
     H.AllPoints = gobjects(0); % Indicate no single handle for all points

     ax_h = struct2cell(H.GroupPoints);
     isPlotted = cellfun(@(x) ~isempty(x), ax_h);
     legend([ax_h{isPlotted}], groupNames{isPlotted}, "Location", "best");
     
else
    % No annotation, plot all points with default settings
    fprintf('Plotting all %d channels with default settings.\n', nChannels);
    h_pts = scatter3(ax, x, y, z, options.MarkerSize, options.DefaultMarkerColor, 'filled');
    H.AllPoints = h_pts;
    % Initialize other handle fields as empty graphics objects
    H.DefaultPoints = gobjects(0);
    H.GroupPoints = struct();
    H.GroupLabels = struct();
end


% --- Final Formatting ---
axis(ax, 'equal');          % Ensure correct aspect ratio
grid(ax, 'on');             % Turn on grid
xlabel(ax, 'X');            % Label axes
ylabel(ax, 'Y');
zlabel(ax, 'Z');
title(ax, '3D Electrode Locations'); % Add title
rotate3d(ax, 'on');         % Enable interactive rotation

% Add lighting if a surface exists and was plotted successfully
if isfield(H, 'Surface') && isgraphics(H.Surface)
    camlight(ax, 'left');
    lighting(ax, 'phong'); % Or 'gouraud'
end

hold(ax, 'off');            % Release the hold on the axes
fprintf('Plotting complete.\n');

end % End of main function topoplot3d

% --- Custom Validation Function for Colors ---
function validatecolor(color)
    if ischar(color) || isstring(color)
        % Check if it's a valid single character color spec
        validChars = 'rgbcmykw';
        if ~ismember(lower(char(color)), validChars)
            error('Invalid color character. Must be one of ''%s''.', validChars);
        end
    elseif isnumeric(color)
        % Check if it's a valid RGB triplet
        if ~isequal(size(color), [1, 3]) || any(color < 0) || any(color > 1)
            error('Invalid RGB triplet. Must be a 1x3 vector with values between 0 and 1.');
        end
    else
        error('Color must be a valid character (e.g., ''k'') or an RGB triplet (e.g., [0 0 1]).');
    end
end
