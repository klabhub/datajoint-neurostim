classdef PupilTracker < handle
    % A class for pupil tracking in Scanbox experiments.
    % Usage:
    % 1.  Create an instance by passing an ns.Experiment table,
    %      the class will locate the associated _eye.mj2 files.
    % 2. Run initialize to manually select the location of the pupil and
    %       the eye. Pupil selection is used as a starting point and to
    %       determine the typical brightness in the pupil. Eye selection is
    %       used to avoid searching for the pupil in impossible locations.
    %      The output of this process is saved to a JSON file next to the
    %      mj2 movie and can be re-read at a later time. 
    % 3. Run track. This will use the initialization parameters and do
    %       automated detection for the full movie. Results are saved to a
    %       .tsv file with a BIDS format json sidecar (both in the same
    %       folder as the .mj2).
    %
    % BK - Mar 2026
    properties (Constant)
        MaxSampleFrames = 50  % Frames used in the interactive initialization routine
        FigureTag = 'sbx.PupilTracker.Figure'  % Unique name of the reused figure
    end

    properties
        Experiment % The ns.Experiment  table
        Files       % A query table of the eye tracking files (ns.File)
        Parameters  % Table of parameters for tracking, determined with the interactive initialization routine.
        Results     % Cell array of tables with tracking results for each video
        Figure      % Handle to the GUI figure
    end

    methods
        function obj = PupilTracker(experiment, pv)
            % PUPIL Construct an instance of the pupil tracking class
            %   obj = sbx.PupilTracker(experiment)
            arguments
                experiment (1,1) ns.Experiment
                pv.filename (1,1) string = 'filename LIKE "%eye.mj2"'
            end
            obj.Experiment = experiment;
            %$ Find sessions with eye tracking movie in the scanbox setup
            obj.Files = obj.Experiment * (ns.File & pv.filename);
            % Pre-populate parameters table
            fileList = fullfile(folder(ns.Experiment & proj(obj.Files)), fetchn(obj.Files, 'filename'));
            numFiles = numel(fileList);
            obj.Parameters = table('Size', [numFiles, 11], ...
                'VariableTypes', {'string', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double'}, ...
                'VariableNames', {'FileName', 'Threshold', 'StartX', 'StartY', 'SeRadius', 'MinArea', 'EyeCenterX', 'EyeCenterY', 'EyeMajorAxis', 'EyeMinorAxis', 'EyeOrientation'});
            obj.Parameters.FileName = fileList;

            %% Check for existing output files
            outputFiles = strrep(obj.Parameters.FileName, '.mj2', '_pupil.tsv');
            fileExists = false(numFiles, 1);
            for i = 1:numFiles
                fileExists(i) = exist(outputFiles(i), "file");
            end

            if any(fileExists)
                fprintf('\nWarning: %d output file(s) already exist:\n', sum(fileExists));
                disp(outputFiles(fileExists));
                response = input('Overwrite, Skip, Read, or Cancel? (o/s/r/c): ', 's');
                switch lower(response)
                    case 'o'
                        % Overwrite: Do nothing, the track method will overwrite
                    case 's'
                        % Skip: Remove existing files from the parameters table
                        obj.Parameters(fileExists, :) = [];
                        fprintf('Skipping %d existing files.\n', sum(fileExists));
                    case 'r'
                        % Read: Load existing TSV files into Results, remove from Parameters
                        obj.Results = cell(numFiles, 1);
                        for i = find(fileExists(:))'
                            obj.Results{i} = readtable(outputFiles(i), 'FileType', 'text', 'Delimiter', '\t');
                            fprintf('  -> Loaded: %s\n', outputFiles(i));
                        end
                        obj.Parameters(fileExists, :) = [];
                        fprintf('Read %d existing file(s). Remaining %d file(s) queued for tracking.\n', ...
                            sum(fileExists), height(obj.Parameters));
                    case 'c'
                        % Cancel: Throw an error and stop
                        error('Execution cancelled by user.');
                    otherwise
                        error('Invalid input. Please choose ''o'', ''s'', ''r'', or ''c''.');
                end
            end
        end

        function disp(obj)
            % DISP - Display a summary of the object.
            fprintf('sbx.PupilTracker with %d file(s) queued for processing:\n', height(obj.Parameters));
            for i = 1:height(obj.Parameters)
                fprintf('  [%d] %s\n', i, obj.Parameters.FileName(i));
            end
        end

        function delete(obj)
            % DELETE - Close the figure when the object is deleted.
            if ~isempty(obj.Figure) && ishandle(obj.Figure)
                close(obj.Figure);
            end
        end

        function initialize(obj, pv)
            % INITIALIZE - Interactively define tracking parameters for multiple videos.
            arguments
                obj
                pv.overwrite (1,1) logical = false
            end

            disp('Starting batch parameter initialization...');
            for i = 1:height(obj.Parameters)
                currentFile = obj.Parameters.FileName(i);
                jsonFileName = strrep(currentFile, '.mj2', '_params.json');
                 % Check whether this exists already
                if exist(jsonFileName, 'file') && ~pv.overwrite
                    response = input(sprintf('Parameter file already exists for %s. Overwrite, Skip, or Read? (o/s/r): ', currentFile), 's');
                    switch lower(response)
                        case 'o'
                            % Proceed with interactive definition
                        case 's'
                            fprintf('Skipping parameter definition for %s.\n', currentFile);
                            continue;
                        case 'r'
                            fprintf('Reading parameters from %s.\n', jsonFileName);
                            try
                                paramStruct = jsondecode(fileread(jsonFileName));
                                obj.Parameters(i, :) = struct2table(paramStruct);
                            catch
                                warning('Could not read or parse %s. Proceeding with interactive definition.', jsonFileName);
                            end
                            continue; % Move to the next file
                        otherwise
                            warning('Invalid input. Proceeding with interactive definition.');
                    end
                end

                fprintf('Processing %d of %d: %s\n', i, height(obj.Parameters), currentFile);
                try
                    v = VideoReader(char(currentFile));
                catch
                    warning('Could not read file. Skipping...');
                    continue;
                end

                % Get random frame indices to sample typical pupil images
                totalFrames = round(v.Duration * v.FrameRate);
                numSample = min(sbx.PupilTracker.MaxSampleFrames, ceil(totalFrames / 20)); % Sample ~5% or max 50 frames
                randomFrameIdx = sort(randperm(totalFrames, numSample));

                % Read random frames and score brightness
                fprintf('Sampling frames for brightness scoring...\n');
                img = zeros(v.Height, v.Width, numSample, sprintf('uint%d', v.BitsPerPixel));
                for ix = 1:numSample
                    try
                        img(:,:,ix) = read(v, randomFrameIdx(ix));
                    catch
                        % The max number of frames is estimated so an error
                        % is possible. Just ignore.
                    end
                end

                brightness = squeeze(median(img, [1 2]));
                [~, sortedIdx] = sort(brightness, 'descend');
                numToAverage = min(10, numSample);
                topFramesIdx = sortedIdx(1:numToAverage);
                % Average the top 10 for visualization
                avgFrameGray = uint8(mean(img(:,:,topFramesIdx), 3));

                % --- Interactive Input ---
                obj.Figure = findobj('Type', 'figure', 'Tag', obj.FigureTag);
                if isempty(obj.Figure)
                    obj.Figure = figure('Name', 'Pupil Initialization', 'Tag', obj.FigureTag, 'units', 'normalized', 'position', [0.25 0.25 0.5 0.5]);
                else
                    figure(obj.Figure); % Bring to front
                    clf(obj.Figure);
                    obj.Figure.Name = 'Pupil Initialization';
                end
                
                imshow(avgFrameGray, 'InitialMagnification', 'fit');
                title(sprintf('Draw an ellipse where the bright pupil appears.\nDouble-click inside to confirm.\nFile: %s', currentFile), 'Interpreter', 'none');
                pupilRoi = drawellipse('Color', 'r', 'FaceAlpha', 0.2);
                wait(pupilRoi);

                title(sprintf('Draw an ellipse for the outer boundary of the eye.\nDouble-click to confirm.'));
                eyeRoi = drawellipse('Color', 'b', 'FaceAlpha', 0.2);
                wait(eyeRoi);

                % --- Extract and Calculate Parameters ---
                obj.Parameters.StartX(i) = pupilRoi.Center(1);
                obj.Parameters.StartY(i) = pupilRoi.Center(2);
                
                pupilMask = createMask(pupilRoi);
                eyeMask = createMask(eyeRoi);
                
                % Calculate threshold from all sampled frames in the masked region        
                eyeNotPupil = img(repmat(eyeMask & ~pupilMask, [1, 1, size(img, 3)]));
                eyeNotPupilLevel = mean(eyeNotPupil, "all");
                pupil = img(repmat(pupilMask, [1, 1, size(img, 3)]));
                pupilLevel = mean(pupil,"all");
                obj.Parameters.Thresholds(i) = (eyeNotPupilLevel + pupilLevel) / 2; % Midpoint between pupil and surrounding eye brightness
        
                     
                % Eye Ellipse parameters
                obj.Parameters.EyeCenterX(i) = eyeRoi.Center(1);
                obj.Parameters.EyeCenterY(i) = eyeRoi.Center(2);
                obj.Parameters.EyeMajorAxis(i) = eyeRoi.SemiAxes(1) * 2;
                obj.Parameters.EyeMinorAxis(i) = eyeRoi.SemiAxes(2) * 2;
                obj.Parameters.EyeOrientation(i) = eyeRoi.RotationAngle;

                % Dynamic Parameter Estimation
                drawnArea = sum(pupilMask(:));
                obj.Parameters.SeRadius(i) = max(2, round(sqrt(drawnArea / pi) * 0.1)); % 10% of radius (min 2)
                obj.Parameters.MinArea(i) = max(15, round(drawnArea * 0.05)); % 5% of drawn area (min 15px)

                % Do not close the figure, just clear it for the next one
                clf(obj.Figure);

                % Write parameters to a JSON file (store only filename, not full path)
                paramStruct = table2struct(obj.Parameters(i, :));
                paramStruct.FileName = char(obj.Parameters.FileName(i));
                [~, paramStruct.FileName, ext] = fileparts(paramStruct.FileName);
                paramStruct.FileName = [paramStruct.FileName ext];
                jsonTxt = jsonencode(paramStruct, "PrettyPrint", true);
                fid = fopen(jsonFileName, 'w');
                if fid == -1
                    warning('Could not write parameters to %s', jsonFileName);
                else
                    fwrite(fid, jsonTxt, 'char');
                    fclose(fid);
                    fprintf('  -> Saved parameters to: %s\n', jsonFileName);
                end
            end
            disp('Finished collecting parameters. Run track() to generate a pupil tracking tsv file.');
        end

        function track(obj, pv)
            % TRACK - Run pupil tracking on all videos.
            arguments
                obj
                pv.minThreshold (1,1) double = 0.1 % Fraction of the initial threshold to allow when adjusting
                pv.visualize (1,1) logical = false  % Set to true to show the frames and the detected regions
                pv.maxFrames (1,1) double = inf  % For quick debugging; process only up to maxFrames
            end

            if isempty(obj.Parameters)
                disp('All files were skipped. Nothing to process.');
                return;
            end

            % Check if parameters have been initialized (Threshold is 0 for all rows
            % when initialize() has not yet been run).
            if all(obj.Parameters.Threshold == 0)
                fprintf('No parameters found. Running initialize() first...\n');
                obj.initialize();
            end

            numVideos = height(obj.Parameters);
            fprintf('Starting ellipse and blob tracking for %d videos...\n', numVideos);

            outputFiles = strrep(obj.Parameters.FileName, '.mj2', '_pupil.tsv');

            if pv.visualize
                obj.Figure = findobj('Type', 'figure', 'Tag', obj.FigureTag);
                if isempty(obj.Figure)
                    obj.Figure = figure('Name', 'Pupil Tracker Viewer', 'Tag', obj.FigureTag, 'Units', 'Normalized', 'Position', [0.25 0.25 0.5 0.5]);
                else
                    figure(obj.Figure); % Bring to front
                    clf(obj.Figure);
                    obj.Figure.Name = 'Pupil Tracker Viewer';
                end
            end

            obj.Results = cell(numVideos, 1);

            for vIdx = 1:numVideos
                videoFile = obj.Parameters.FileName(vIdx);

                % Load Dynamic Parameters
                pupilThreshold = obj.Parameters.Threshold(vIdx);
                se = strel('disk', obj.Parameters.SeRadius(vIdx));
                minPupilArea = obj.Parameters.MinArea(vIdx);

                % Eye boundary parameters
                eyeCenterX = obj.Parameters.EyeCenterX(vIdx);
                eyeCenterY = obj.Parameters.EyeCenterY(vIdx);
                eyeMajorAxis = obj.Parameters.EyeMajorAxis(vIdx);
                eyeMinorAxis = obj.Parameters.EyeMinorAxis(vIdx);
                eyeOrientation = obj.Parameters.EyeOrientation(vIdx);

                fprintf('\n[%d/%d] Processing: %s\n \t\t -> %s\n', vIdx, numVideos, videoFile, outputFiles(vIdx));

                try
                    v = VideoReader(char(videoFile));
                catch
                    warning('Could not read %s. Skipping...', videoFile);
                    continue;
                end

                estimatedFrames = ceil(v.Duration * v.FrameRate);

                % Core Arrays
                Frames = (1:estimatedFrames)';
                Centroid_X = NaN(estimatedFrames, 1);
                Centroid_Y = NaN(estimatedFrames, 1);
                PupilArea = NaN(estimatedFrames, 1);
                ThresholdUsed = NaN(estimatedFrames, 1);

                % Ellipse-specific Arrays
                MajorAxis = NaN(estimatedFrames, 1);
                MinorAxis = NaN(estimatedFrames, 1);
                Eccentricity = NaN(estimatedFrames, 1);
                Orientation = NaN(estimatedFrames, 1);
                EM = NaN(estimatedFrames, 1);
                FitQuality = NaN(estimatedFrames, 1);
                IntensityRatio = NaN(estimatedFrames, 1);

                % Blob-specific Arrays
                BoundingBox_X = NaN(estimatedFrames, 1);
                BoundingBox_Y = NaN(estimatedFrames, 1);
                BoundingBox_Width = NaN(estimatedFrames, 1);
                BoundingBox_Height = NaN(estimatedFrames, 1);

                [Xgrid, Ygrid] = meshgrid(1:v.Width, 1:v.Height);

                % Create a static elliptical mask for the eye
                a = eyeMajorAxis / 2;
                b = eyeMinorAxis / 2;
                theta = -deg2rad(eyeOrientation); % Convert to radians and invert for calculation

                % Rotated ellipse equation
                eyeMask = (((Xgrid - eyeCenterX) * cos(theta) + (Ygrid - eyeCenterY) * sin(theta)).^2 / a^2 + ...
                    ((Xgrid - eyeCenterX) * sin(theta) - (Ygrid - eyeCenterY) * cos(theta)).^2 / b^2) <= 1;

                frameCount = 0;

                while hasFrame(v) && frameCount <= pv.maxFrames
                    frameCount = frameCount + 1;
                    frame = readFrame(v);
                    if size(frame, 3) == 3, frame = rgb2gray(frame); end

                    [level, em] = graythresh(frame(eyeMask));
                    EM(frameCount) = em;
                    % Quality check on graythresh (effectiveness metric) and threshold limit
                    if em > 0.5 % Check if histogram is reasonably bimodal
                        currentThreshold = level * max(frame(:));
                        % Enforce minimum threshold to be at least a fraction of the original
                        currentThreshold = max(currentThreshold, pv.minThreshold* pupilThreshold);
                    else
                        % graythresh failed, likely no pupil, use default
                        currentThreshold = pupilThreshold;
                    end

                    if pv.visualize  && ishandle(obj.Figure)
                        imshow(frame, 'Parent', gca(obj.Figure),'initialMagnification','fit'); hold(gca(obj.Figure), 'on');
                    end
                    % Image Processing: Search within the defined eye boundary
                    binaryImage = (frame > currentThreshold) & eyeMask;
                    cleanImage = bwareaopen(imopen(imfill(binaryImage, 'holes'), se), minPupilArea);
                    ThresholdUsed(frameCount) = currentThreshold;

                    % Calculate intensity ratio quality metric
                    if any(cleanImage(:))
                        insideMask = cleanImage;
                        outsideMask = eyeMask & ~insideMask;
                        meanInside = mean(frame(insideMask));
                        meanOutside = mean(frame(outsideMask));
                        if meanInside > 0
                            IntensityRatio(frameCount) =  meanInside/meanOutside ; % Higher ratio suggests better contrast
                        end
                    end

                    stats_ellipse = regionprops(cleanImage, 'Area', 'Centroid', 'MajorAxisLength', 'MinorAxisLength', 'Eccentricity', 'Orientation');
                    stats_blob = regionprops(cleanImage, 'Area', 'Centroid', 'BoundingBox');

                    if ~isempty(stats_ellipse)
                        [~, largestIdx] = max([stats_ellipse.Area]);
                        bestEllipse = stats_ellipse(largestIdx);

                        Centroid_X(frameCount) = bestEllipse.Centroid(1);
                        Centroid_Y(frameCount) = bestEllipse.Centroid(2);
                        PupilArea(frameCount) = bestEllipse.Area;

                        MajorAxis(frameCount) = bestEllipse.MajorAxisLength;
                        MinorAxis(frameCount) = bestEllipse.MinorAxisLength;
                        Eccentricity(frameCount) = bestEllipse.Eccentricity;
                        Orientation(frameCount) = bestEllipse.Orientation;

                        ellipseArea = pi * (bestEllipse.MajorAxisLength / 2) * (bestEllipse.MinorAxisLength / 2);
                        if ellipseArea > 0
                            FitQuality(frameCount) = bestEllipse.Area / ellipseArea;
                        end
                    end

                    if ~isempty(stats_blob)
                        [~, largestIdx] = max([stats_blob.Area]);
                        bestBlob = stats_blob(largestIdx);

                        BoundingBox_X(frameCount) = bestBlob.BoundingBox(1);
                        BoundingBox_Y(frameCount) = bestBlob.BoundingBox(2);
                        BoundingBox_Width(frameCount) = bestBlob.BoundingBox(3);
                        BoundingBox_Height(frameCount) = bestBlob.BoundingBox(4);
                    end

                    if pv.visualize && ishandle(obj.Figure)
                        phi_eye = linspace(0, 2*pi, 50);
                        R_eye = [cos(-theta) sin(-theta); -sin(-theta) cos(-theta)];
                        xy_eye = R_eye * [a*cos(phi_eye); b*sin(phi_eye)];
                        plot(obj.Figure.CurrentAxes, xy_eye(1,:) + eyeCenterX, xy_eye(2,:) + eyeCenterY, 'b--', 'LineWidth', 1);

                        if ~isempty(stats_ellipse)
                            plot(obj.Figure.CurrentAxes, bestEllipse.Centroid(1), bestEllipse.Centroid(2), 'g+', 'MarkerSize', 8, 'LineWidth', 2);
                            phi_pupil = linspace(0, 2*pi, 50);
                            a_pupil = bestEllipse.MajorAxisLength / 2;
                            b_pupil = bestEllipse.MinorAxisLength / 2;
                            theta_pupil = pi * bestEllipse.Orientation / 180;
                            R_pupil = [cos(theta_pupil) sin(theta_pupil); -sin(theta_pupil) cos(theta_pupil)];
                            xy_pupil = R_pupil * [a_pupil*cos(phi_pupil); b_pupil*sin(phi_pupil)];
                            plot(obj.Figure.CurrentAxes, xy_pupil(1,:) + bestEllipse.Centroid(1), xy_pupil(2,:) + bestEllipse.Centroid(2), 'r-', 'LineWidth', 2);
                        end

                        if ~isempty(stats_blob)
                            rectangle(obj.Figure.CurrentAxes, 'Position', bestBlob.BoundingBox, 'EdgeColor', 'g', 'LineWidth', 2);
                        end

                        if isempty(stats_blob) && isempty(stats_ellipse)
                            imshow(frame, 'Parent', gca(obj.Figure), 'initialMagnification', 'fit');
                            title(gca(obj.Figure), sprintf('Frame %d | BLINK', frameCount), 'Color', 'r');
                        else
                            title(gca(obj.Figure), sprintf('Frame %d | Area: %d | EM: %.2f | Fit: %.2f | IR: %.2f', frameCount, bestEllipse.Area, em, FitQuality(frameCount), IntensityRatio(frameCount)), 'Interpreter', 'none');
                        end
                        hold(gca(obj.Figure), 'off');
                        drawnow limitrate;
                    end
                end
                
                % Store results as a member variable
                idx = 1:frameCount;
                obj.Results{vIdx} = table(Frames(idx), Centroid_X(idx), Centroid_Y(idx), PupilArea(idx), ...
                    MajorAxis(idx), MinorAxis(idx), Eccentricity(idx), Orientation(idx), ...
                    BoundingBox_X(idx), BoundingBox_Y(idx), BoundingBox_Width(idx), BoundingBox_Height(idx), ...
                    ThresholdUsed(idx), EM(idx), FitQuality(idx), IntensityRatio(idx));
                obj.Results{vIdx}.Properties.VariableNames = {'Frame', 'Centroid_X', 'Centroid_Y', 'Area', ...
                    'MajorAxis', 'MinorAxis', 'Eccentricity', 'Orientation', ...
                    'BBox_X', 'BBox_Y', 'BBox_Width', 'BBox_Height', 'Threshold', 'EM', 'FitQuality', 'IntensityRatio'};

                writetable(obj.Results{vIdx}, outputFiles(vIdx), 'FileType', 'text', 'Delimiter', '\t');
                fprintf('    -> Saved data to: %s\n', outputFiles(vIdx));

                jsonOutputName = strrep(outputFiles(vIdx), '.tsv', '.json');
                sbx.PupilTracker.writeBidsJson(jsonOutputName);
                fprintf('    -> Wrote BIDS sidecar to: %s\n', jsonOutputName);
            end

            if pv.visualize && ishandle(obj.Figure), clf(obj.Figure); end
            disp('Batch processing complete!');
        end

        function read(obj)
            % READ - Load existing _pupil.tsv files into obj.Results.
            numVideos = height(obj.Parameters);
            outputFiles = strrep(obj.Parameters.FileName, '.mj2', '_pupil.tsv');
            obj.Results = cell(numVideos, 1);
            for vIdx = 1:numVideos
                tsvFile = outputFiles(vIdx);
                if exist(tsvFile, 'file')
                    obj.Results{vIdx} = readtable(tsvFile, 'FileType', 'text', 'Delimiter', '\t');
                    fprintf('  -> Loaded: %s\n', tsvFile);
                else
                    fprintf('  -> Not found, skipping: %s\n', tsvFile);
                end
            end
            fprintf('Done. %d/%d result(s) loaded.\n', sum(~cellfun(@isempty, obj.Results)), numVideos);
        end

        function plot(obj, columns, pv)
            % PLOT - Plot one or more result columns for each tracked video.
            %   obj.plot("Area")
            %   obj.plot(["Area","EM","FitQuality"])            
            arguments
                obj
                columns (1,:) string {mustBeMember(columns, ["Frame","Centroid_X","Centroid_Y","Area","MajorAxis","MinorAxis", ...
                                                    "Eccentricity","Orientation","BBox_X","BBox_Y","BBox_Width","BBox_Height", ...
                                                    "Threshold","EM","FitQuality","IntensityRatio"])} = "Area" 
                pv.Parent = []
            end
            validIdx = find(~cellfun(@isempty, obj.Results));
            if isempty(validIdx)
                disp('No results to plot. Run track() or read() first.');
                return;
            end
            if isempty(pv.Parent)
                fig = figure('Name', 'PupilTracker Results');
            else
                fig = pv.Parent;
            end
            tiledlayout(fig, numel(validIdx), 1, 'TileSpacing', 'compact', 'Padding', 'compact');
            for i = 1:numel(validIdx)
                vIdx = validIdx(i);
                t = obj.Results{vIdx};
                ax = nexttile;
                hold(ax, 'on');
                for c = 1:numel(columns)
                    col = columns(c);
                    if ismember(col, t.Properties.VariableNames)
                        plot(ax, t.Frame, t.(col), 'DisplayName', col);
                    else
                        warning('Column "%s" not found in Results{%d}.', col, vIdx);
                    end
                end
                hold(ax, 'off');
                if vIdx <= height(obj.Parameters)
                    [~, fname] = fileparts(obj.Parameters.FileName(vIdx));
                else
                    fname = sprintf('Result %d', vIdx);
                end
                title(ax, fname, 'Interpreter', 'none');
                xlabel(ax, 'Frame');
                if numel(columns) > 1
                    legend(ax, 'Location', 'best');
                else
                    ylabel(ax, columns(1));
                end
            end
        end
    end

    methods (Access = private, Static)
        function writeBidsJson(fileName)
            bidsInfo.Frame = struct('Description', 'Frame number from the video file.', 'Units', 'integer');
            bidsInfo.Centroid_X = struct('Description', 'X-coordinate of the pupil centroid.', 'Units', 'pixels');
            bidsInfo.Centroid_Y = struct('Description', 'Y-coordinate of the pupil centroid.', 'Units', 'pixels');
            bidsInfo.Area = struct('Description', 'Area of the detected pupil.', 'Units', 'pixels^2');
            bidsInfo.MajorAxis = struct('Description', 'Length of the major axis of the ellipse fitted to the pupil.', 'Units', 'pixels');
            bidsInfo.MinorAxis = struct('Description', 'Length of the minor axis of the ellipse fitted to the pupil.', 'Units', 'pixels');
            bidsInfo.Eccentricity = struct('Description', 'Eccentricity of the ellipse fitted to the pupil (0=circle, <1=ellipse).', 'Units', 'ratio');
            bidsInfo.Orientation = struct('Description', 'Orientation of the ellipse fitted to the pupil.', 'Units', 'degrees');
            bidsInfo.BBox_X = struct('Description', 'X-coordinate of the upper-left corner of the bounding box around the pupil.', 'Units', 'pixels');
            bidsInfo.BBox_Y = struct('Description', 'Y-coordinate of the upper-left corner of the bounding box around the pupil.', 'Units', 'pixels');
            bidsInfo.BBox_Width = struct('Description', 'Width of the bounding box around the pupil.', 'Units', 'pixels');
            bidsInfo.BBox_Height = struct('Description', 'Height of the bounding box around the pupil.', 'Units', 'pixels');
            bidsInfo.Threshold = struct('Description', 'The grayscale threshold value used for pupil detection in this frame.', 'Units', 'grayscale_value');
            bidsInfo.EM = struct('Description', 'Effectiveness metric of the threshold, indicating bimodality of the histogram.', 'Units', 'ratio');
            bidsInfo.FitQuality = struct('Description', 'Goodness of fit for the ellipse, calculated as the ratio of pixel area to the geometric area of the fitted ellipse. A value near 1 indicates a good fit.', 'Units', 'ratio');
            bidsInfo.IntensityRatio = struct('Description', 'Ratio of the mean intensity outside the pupil to inside the pupil, within the eye mask. Higher values suggest better contrast.', 'Units', 'ratio');

            jsonTxt = jsonencode(bidsInfo, "PrettyPrint", true);
            fid = fopen(fileName, 'w');
            if fid == -1, error('Cannot create JSON file: %s', fileName); end
            fwrite(fid, jsonTxt, 'char');
            fclose(fid);
        end
    end
end
