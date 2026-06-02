classdef PupilTrackerApp < matlab.apps.AppBase
    % sbx.PupilTrackerApp  App Designer-style interface for pupil setup and tracking.

    properties (Access = public)
        UIFigure
        MainGrid
        SidebarGrid
        WorkspaceTabs
        FileListLabel
        FileListBox
        InitializeButton
        NewAverageButton
        ConfirmEyeButton
        SaveParametersButton
        RestartButton
        SkipButton
        TrackCurrentButton
        TrackCurrentPreviewButton
        TrackAllButton
        TrackAllPreviewButton
        StatusTextArea
        EyeTab
        PupilTab
        PreviewTab
        SummaryTab
        EyeAxes
        PupilGrid
        PreviewAxes
        PreviewSlider
        PreviewLabel
        SummaryTable
    end

    properties (Access = private)
        Tracker
        FileKeys = strings(0,1)
        CurrentFile = ""
        CurrentVideoReader = []
        CurrentVideoFile = ""
        AverageFrame = []
        EyeRoi = []
        PupilRois = gobjects(0)
        PupilAxes = gobjects(0)
        SampleFramePool = []
        PupilFramePool = []
        QuantPositions = []
        EyeMask = []
        EyeCenter = []
        EyeSemiAxes = []
        EyeAngle = 0
        ReviewVideoFile = ""
        ReviewResultTable = table()
        ReviewParameters = table()
        ReviewLastFrame = 0
        ReviewEyeOutline = zeros(2, 0)
        IgnorePreviewSlider = false
    end

    methods (Access = public)
        function app = PupilTrackerApp(tracker, mode)
            arguments
                tracker (1,1) sbx.PupilTracker
                mode (1,1) string {mustBeMember(mode, ["initialize", "track"])} = "initialize"
            end

            app.Tracker = tracker;
            app.createComponents();
            registerApp(app, app.UIFigure);
            app.populateFileList();
            app.switchMode(mode);
            if ~isempty(app.FileKeys)
                app.selectFile(app.FileKeys(1));
            else
                app.appendStatus('No files are currently queued in this PupilTracker object.');
            end
            if nargout == 0
                clear app
            end
        end

        function delete(app)
            if ~isempty(app.UIFigure) && isvalid(app.UIFigure)
                delete(app.UIFigure);
            end
        end

        function switchMode(app, mode)
            if mode == "track"
                app.WorkspaceTabs.SelectedTab = app.PreviewTab;
                app.appendStatus('Track mode selected. Use Track Current or Track All.');
            else
                app.WorkspaceTabs.SelectedTab = app.EyeTab;
                app.appendStatus('Initialize mode selected. Load an average frame, draw the eye ROI, then save pupil parameters.');
            end
        end
    end

    methods (Access = private)
        function createComponents(app)
            app.UIFigure = uifigure('Name', 'Pupil Tracker App', 'Position', [100 100 1500 920]);
            app.UIFigure.CloseRequestFcn = @(~,~) delete(app);

            app.MainGrid = uigridlayout(app.UIFigure, [1 2]);
            app.MainGrid.ColumnWidth = {280, '1x'};
            app.MainGrid.RowHeight = {'1x'};
            app.MainGrid.Padding = [8 8 8 8];

            app.SidebarGrid = uigridlayout(app.MainGrid, [14 1]);
            app.SidebarGrid.Layout.Row = 1;
            app.SidebarGrid.Layout.Column = 1;
            app.SidebarGrid.RowHeight = {22, 170, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, '1x', 180};
            app.SidebarGrid.ColumnWidth = {'1x'};

            app.FileListLabel = uilabel(app.SidebarGrid, 'Text', 'Queued Files');
            app.FileListLabel.Layout.Row = 1;

            app.FileListBox = uilistbox(app.SidebarGrid, 'ValueChangedFcn', @(~,~) app.onFileChanged());
            app.FileListBox.Layout.Row = 2;

            app.InitializeButton = uibutton(app.SidebarGrid, 'push', 'Text', 'Load Eye Average', ...
                'ButtonPushedFcn', @(~,~) app.loadEyeAverage());
            app.InitializeButton.Layout.Row = 3;

            app.NewAverageButton = uibutton(app.SidebarGrid, 'push', 'Text', 'New Average', ...
                'ButtonPushedFcn', @(~,~) app.loadEyeAverage());
            app.NewAverageButton.Layout.Row = 4;

            app.ConfirmEyeButton = uibutton(app.SidebarGrid, 'push', 'Text', 'Confirm Eye / Load Pupils', ...
                'ButtonPushedFcn', @(~,~) app.confirmEye());
            app.ConfirmEyeButton.Layout.Row = 5;

            app.SaveParametersButton = uibutton(app.SidebarGrid, 'push', 'Text', 'Save Parameters', ...
                'ButtonPushedFcn', @(~,~) app.saveParameters());
            app.SaveParametersButton.Layout.Row = 6;

            app.RestartButton = uibutton(app.SidebarGrid, 'push', 'Text', 'Restart Current File', ...
                'ButtonPushedFcn', @(~,~) app.restartCurrentFile());
            app.RestartButton.Layout.Row = 7;

            app.SkipButton = uibutton(app.SidebarGrid, 'push', 'Text', 'Skip Current File', ...
                'ButtonPushedFcn', @(~,~) app.skipCurrentFile());
            app.SkipButton.Layout.Row = 8;

            app.TrackCurrentPreviewButton = uibutton(app.SidebarGrid, 'push', 'Text', 'Track Current Preview', ...
                'ButtonPushedFcn', @(~,~) app.trackFiles(false, true));
            app.TrackCurrentPreviewButton.Layout.Row = 9;

            app.TrackCurrentButton = uibutton(app.SidebarGrid, 'push', 'Text', 'Track Current', ...
                'ButtonPushedFcn', @(~,~) app.trackFiles(false, false));
            app.TrackCurrentButton.Layout.Row = 10;

            app.TrackAllPreviewButton = uibutton(app.SidebarGrid, 'push', 'Text', 'Track All Preview', ...
                'ButtonPushedFcn', @(~,~) app.trackFiles(true, true));
            app.TrackAllPreviewButton.Layout.Row = 11;

            app.TrackAllButton = uibutton(app.SidebarGrid, 'push', 'Text', 'Track All', ...
                'ButtonPushedFcn', @(~,~) app.trackFiles(true, false));
            app.TrackAllButton.Layout.Row = 12;

            app.StatusTextArea = uitextarea(app.SidebarGrid, 'Editable', 'off');
            app.StatusTextArea.Layout.Row = 14;
            app.StatusTextArea.Value = {'Ready.'};

            app.WorkspaceTabs = uitabgroup(app.MainGrid);
            app.WorkspaceTabs.Layout.Row = 1;
            app.WorkspaceTabs.Layout.Column = 2;

            app.EyeTab = uitab(app.WorkspaceTabs, 'Title', 'Eye ROI');
            app.EyeAxes = uiaxes(app.EyeTab, 'Position', [10 10 1180 820]);
            app.configureImageAxes(app.EyeAxes);

            app.PupilTab = uitab(app.WorkspaceTabs, 'Title', 'Pupil ROIs');
            app.PupilGrid = uigridlayout(app.PupilTab, [2 3]);
            app.PupilGrid.RowHeight = {'1x', '1x'};
            app.PupilGrid.ColumnWidth = {'1x', '1x', '1x'};
            app.PupilGrid.Padding = [8 8 8 8];
            app.PupilGrid.RowSpacing = 8;
            app.PupilGrid.ColumnSpacing = 8;
            app.PupilAxes = gobjects(app.Tracker.NrPupilImages, 1);
            for idx = 1:app.Tracker.NrPupilImages
                app.PupilAxes(idx) = uiaxes(app.PupilGrid);
                app.configureImageAxes(app.PupilAxes(idx));
            end

            app.PreviewTab = uitab(app.WorkspaceTabs, 'Title', 'Tracking Preview');
            previewGrid = uigridlayout(app.PreviewTab, [3 1]);
            previewGrid.RowHeight = {'1x', 40, 30};
            previewGrid.ColumnWidth = {'1x'};
            previewGrid.Padding = [8 8 8 8];
            previewGrid.RowSpacing = 8;

            app.PreviewAxes = uiaxes(previewGrid);
            app.PreviewAxes.Layout.Row = 1;
            app.configureImageAxes(app.PreviewAxes);

            app.PreviewSlider = uislider(previewGrid, ...
                'Limits', [1 2], 'Value', 1, 'Enable', 'off', ...
                'ValueChangedFcn', @(src,~) app.onPreviewSliderChanged(src.Value), ...
                'ValueChangingFcn', @(~,evt) app.onPreviewSliderChanged(evt.Value));
            app.PreviewSlider.Layout.Row = 2;

            app.PreviewLabel = uilabel(previewGrid, 'Text', '0 / 0');
            app.PreviewLabel.Layout.Row = 3;

            app.SummaryTab = uitab(app.WorkspaceTabs, 'Title', 'Summary');
            app.SummaryTable = uitable(app.SummaryTab, 'Position', [10 10 1180 820]);
            app.SummaryTable.ColumnEditable = false;
        end

        function configureImageAxes(~, ax)
            ax.Toolbar.Visible = 'off';
            ax.XTick = [];
            ax.YTick = [];
            ax.Box = 'on';
            ax.YDir = 'reverse';
            ax.DataAspectRatio = [1 1 1];
        end

        function populateFileList(app)
            app.FileKeys = string(keys(app.Tracker.Parameters))';
            if isempty(app.FileKeys)
                app.FileListBox.Items = {};
                app.FileListBox.ItemsData = {};
                return;
            end

            labels = cell(numel(app.FileKeys), 1);
            itemsData = cell(numel(app.FileKeys), 1);
            for idx = 1:numel(app.FileKeys)
                fileKey = app.FileKeys(idx);
                labels{idx} = app.makeFileLabel(fileKey);
                itemsData{idx} = char(fileKey);
            end
            app.FileListBox.Items = labels;
            app.FileListBox.ItemsData = itemsData;
            if isempty(app.FileListBox.Value) || ~any(strcmp(app.FileListBox.Value, itemsData))
                app.FileListBox.Value = itemsData{1};
            end
        end

        function label = makeFileLabel(app, fileKey)
            [~, name, ext] = fileparts(char(fileKey));
            hasParams = isKey(app.Tracker.Parameters, char(fileKey)) && app.Tracker.Parameters(char(fileKey)).Threshold ~= 0;
            hasResults = isKey(app.Tracker.Results, char(fileKey)) && ~isempty(app.Tracker.Results(char(fileKey)));
            if hasParams
                pFlag = 'P';
            else
                pFlag = '-';
            end
            if hasResults
                rFlag = 'R';
            else
                rFlag = '-';
            end
            label = sprintf('[%s%s] %s%s', pFlag, rFlag, name, ext);
        end

        function onFileChanged(app)
            if isempty(app.FileListBox.ItemsData)
                return;
            end
            app.selectFile(string(app.FileListBox.Value));
        end

        function selectFile(app, fileKey)
            if strlength(fileKey) == 0
                return;
            end
            app.CurrentFile = fileKey;
            app.resetTransientState();
            app.updateSummary();
            app.appendStatus("Selected " + fileKey);
            if isKey(app.Tracker.Results, char(fileKey)) && ~isempty(app.Tracker.Results(char(fileKey)))
                app.prepareReviewData(fileKey, app.Tracker.Results(char(fileKey)));
            else
                app.PreviewSlider.Enable = 'off';
                app.PreviewLabel.Text = '0 / 0';
                hold(app.PreviewAxes, 'off');
                cla(app.PreviewAxes);
            end
        end

        function resetTransientState(app)
            if ~isempty(app.EyeRoi) && isvalid(app.EyeRoi)
                delete(app.EyeRoi);
            end
            app.EyeRoi = [];
            if ~isempty(app.PupilRois)
                for idx = 1:numel(app.PupilRois)
                    if isgraphics(app.PupilRois(idx))
                        delete(app.PupilRois(idx));
                    end
                end
            end
            app.PupilRois = gobjects(0);
            app.SampleFramePool = [];
            app.PupilFramePool = [];
            app.QuantPositions = [];
            app.EyeMask = [];
            app.EyeCenter = [];
            app.EyeSemiAxes = [];
            app.EyeAngle = 0;
            hold(app.EyeAxes, 'off');
            cla(app.EyeAxes);
            for idx = 1:numel(app.PupilAxes)
                hold(app.PupilAxes(idx), 'off');
                cla(app.PupilAxes(idx));
            end
        end

        function loadEyeAverage(app)
            if ~app.ensureVideoLoaded()
                return;
            end
            app.SampleFramePool = app.getFramePool(app.CurrentVideoReader, app.Tracker.MaxSampleFrames);
            brightnessPerFrame = squeeze(mean(app.SampleFramePool, [1 2]));
            keep = brightnessPerFrame > prctile(brightnessPerFrame, 50);
            if ~any(keep)
                keep = true(size(brightnessPerFrame));
            end
            app.AverageFrame = uint8(mean(double(app.SampleFramePool(:, :, keep)), 3));
            hold(app.EyeAxes, 'off');
            imshow(app.AverageFrame, 'Parent', app.EyeAxes, 'InitialMagnification', 'fit');
            title(app.EyeAxes, 'Draw eye boundary, then click Confirm Eye / Load Pupils.');
            if ~isempty(app.EyeRoi) && isvalid(app.EyeRoi)
                delete(app.EyeRoi);
            end
            app.EyeRoi = drawellipse(app.EyeAxes, 'Color', 'b', 'FaceAlpha', 0.15);
            app.WorkspaceTabs.SelectedTab = app.EyeTab;
            app.appendStatus('Loaded a new average eye frame.');
        end

        function confirmEye(app)
            if isempty(app.EyeRoi) || ~isvalid(app.EyeRoi)
                uialert(app.UIFigure, 'Draw an eye ellipse first.', 'Eye ROI Required');
                return;
            end
            if isempty(app.SampleFramePool)
                uialert(app.UIFigure, 'Load an average eye frame first.', 'Sample Frames Required');
                return;
            end

            app.EyeCenter = app.EyeRoi.Center;
            app.EyeSemiAxes = app.EyeRoi.SemiAxes;
            app.EyeAngle = app.EyeRoi.RotationAngle;

            [Xg, Yg] = meshgrid(1:app.CurrentVideoReader.Width, 1:app.CurrentVideoReader.Height);
            aEye = app.EyeSemiAxes(1);
            bEye = app.EyeSemiAxes(2);
            thetaEye = -deg2rad(app.EyeAngle);
            app.EyeMask = (((Xg - app.EyeCenter(1)) * cos(thetaEye) + (Yg - app.EyeCenter(2)) * sin(thetaEye)).^2 / aEye^2 + ...
                ((Xg - app.EyeCenter(1)) * sin(thetaEye) - (Yg - app.EyeCenter(2)) * cos(thetaEye)).^2 / bEye^2) <= 1;

            eyePixelsRef = double(app.AverageFrame(app.EyeMask));
            eyeFloorRef = prctile(eyePixelsRef, 5);
            nPool = size(app.SampleFramePool, 3);
            eyeBrightness = zeros(nPool, 1);
            keepFrame = true(nPool, 1);
            for frameIdx = 1:nPool
                frameSlice = double(app.SampleFramePool(:, :, frameIdx));
                eyeBrightness(frameIdx) = mean(frameSlice(app.EyeMask));
                if eyeBrightness(frameIdx) < eyeFloorRef
                    keepFrame(frameIdx) = false;
                end
            end
            keptIdx = find(keepFrame);
            if numel(keptIdx) < app.Tracker.NrPupilImages
                keptIdx = 1:nPool;
            end
            [~, rankAmongKept] = sort(eyeBrightness(keptIdx));
            rankOrder = keptIdx(rankAmongKept);
            nKept = numel(rankOrder);
            app.QuantPositions = round(linspace(1, nKept, app.Tracker.NrPupilImages));
            tileIdx = rankOrder(app.QuantPositions);
            app.PupilFramePool = app.SampleFramePool(:, :, tileIdx);

            phiEye = linspace(0, 2 * pi, 100);
            alpha = deg2rad(app.EyeAngle);
            rotationEye = [cos(alpha) sin(alpha); -sin(alpha) cos(alpha)];
            xyEye = rotationEye * [aEye * cos(phiEye); bEye * sin(phiEye)];

            if ~isempty(app.PupilRois)
                for idx = 1:numel(app.PupilRois)
                    if isgraphics(app.PupilRois(idx))
                        delete(app.PupilRois(idx));
                    end
                end
            end
            app.PupilRois = gobjects(app.Tracker.NrPupilImages, 1);
            for idx = 1:app.Tracker.NrPupilImages
                ax = app.PupilAxes(idx);
                hold(ax, 'off');
                imshow(app.PupilFramePool(:, :, idx), 'Parent', ax, 'InitialMagnification', 'fit');
                hold(ax, 'on');
                plot(ax, xyEye(1, :) + app.EyeCenter(1), xyEye(2, :) + app.EyeCenter(2), 'b-', 'LineWidth', 1.5);
                title(ax, sprintf('Rank %d/%d', app.QuantPositions(idx), nKept));
                app.PupilRois(idx) = drawellipse(ax, 'Color', 'r', 'FaceAlpha', 0.15);
                hold(ax, 'off');
            end
            app.WorkspaceTabs.SelectedTab = app.PupilTab;
            app.appendStatus('Eye ROI confirmed. Draw or adjust the red pupil ellipses, then click Save Parameters.');
        end

        function saveParameters(app)
            if strlength(app.CurrentFile) == 0 || isempty(app.PupilFramePool) || isempty(app.EyeMask)
                uialert(app.UIFigure, 'Complete eye and pupil selection first.', 'Parameters Not Ready');
                return;
            end
            allPupilPixels = [];
            allEyePixels = [];
            validCentroids = zeros(0, 2);
            validAreas = zeros(0, 1);
            for idx = 1:numel(app.PupilRois)
                if ~isgraphics(app.PupilRois(idx))
                    continue;
                end
                cx = app.PupilRois(idx).Center(1);
                cy = app.PupilRois(idx).Center(2);
                row = round(cy);
                col = round(cx);
                if row < 1 || row > app.CurrentVideoReader.Height || col < 1 || col > app.CurrentVideoReader.Width
                    continue;
                end
                if ~app.EyeMask(row, col)
                    continue;
                end
                mask = createMask(app.PupilRois(idx));
                frame = app.PupilFramePool(:, :, idx);
                thisEyePixels = double(frame(app.EyeMask & ~mask));
                thisPupilPixels = double(frame(mask));
                allPupilPixels = [allPupilPixels; thisPupilPixels]; %#ok<AGROW>
                allEyePixels = [allEyePixels; thisEyePixels]; %#ok<AGROW>
                validCentroids(end+1, :) = [cx, cy]; %#ok<AGROW>
                validAreas(end+1) = pi * prod(app.PupilRois(idx).SemiAxes); %#ok<AGROW>
            end

            if isempty(allPupilPixels)
                uialert(app.UIFigure, 'No valid pupil ellipses are inside the eye ROI.', 'No Valid Pupil ROIs');
                return;
            end

            app.Tracker.Parameters(char(app.CurrentFile)) = table( ...
                (median(allPupilPixels) + median(allEyePixels)) / 2, ...
                median(allEyePixels), ...
                median(validCentroids(:, 1)), median(validCentroids(:, 2)), ...
                median(validAreas), ...
                app.EyeCenter(1), app.EyeCenter(2), ...
                app.EyeSemiAxes(1) * 2, app.EyeSemiAxes(2) * 2, app.EyeAngle, ...
                'VariableNames', {'Threshold', 'Floor', 'StartX', 'StartY', 'Area', 'EyeX', 'EyeY', 'EyeMajor', 'EyeMinor', 'EyeOrientation'});

            [folderName, fileName, ext] = fileparts(char(app.CurrentFile));
            paramStruct = table2struct(app.Tracker.Parameters(char(app.CurrentFile)));
            paramStruct.FileName = [fileName ext];
            jsonText = jsonencode(paramStruct, 'PrettyPrint', true);
            jsonFile = fullfile(folderName, strrep([fileName ext], '.mj2', '_pupil.json'));
            fileId = fopen(jsonFile, 'w');
            if fileId == -1
                uialert(app.UIFigure, ['Could not write ' jsonFile], 'Write Error');
                return;
            end
            fwrite(fileId, jsonText, 'char');
            fclose(fileId);

            app.populateFileList();
            app.updateSummary();
            app.appendStatus("Saved parameters to " + jsonFile);
        end

        function restartCurrentFile(app)
            if strlength(app.CurrentFile) == 0
                return;
            end
            app.resetTransientState();
            app.loadEyeAverage();
            app.appendStatus('Restarted current file initialization.');
        end

        function skipCurrentFile(app)
            if strlength(app.CurrentFile) == 0
                return;
            end
            current = char(app.CurrentFile);
            if isKey(app.Tracker.Parameters, current)
                remove(app.Tracker.Parameters, current);
            end
            if isKey(app.Tracker.Results, current)
                remove(app.Tracker.Results, current);
            end
            app.appendStatus("Skipped " + app.CurrentFile);
            app.populateFileList();
            if ~isempty(app.FileKeys)
                app.selectFile(app.FileKeys(1));
            else
                app.CurrentFile = "";
                app.resetTransientState();
                app.SummaryTable.Data = table();
            end
        end

        function trackFiles(app, allFiles, preview)
            if isempty(app.FileKeys)
                uialert(app.UIFigure, 'There are no queued files to track.', 'Nothing To Track');
                return;
            end
            if allFiles
                fileSet = app.FileKeys;
            else
                fileSet = app.CurrentFile;
            end
            fileSet = fileSet(strlength(fileSet) > 0);
            if isempty(fileSet)
                return;
            end
            app.WorkspaceTabs.SelectedTab = app.PreviewTab;
            for idx = 1:numel(fileSet)
                fileKey = fileSet(idx);
                if ~isKey(app.Tracker.Parameters, char(fileKey)) || app.Tracker.Parameters(char(fileKey)).Threshold == 0
                    app.appendStatus("Skipping " + fileKey + ' because parameters are not initialized.');
                    continue;
                end
                app.runTracking(fileKey, preview);
                app.populateFileList();
            end
            app.updateSummary();
        end

        function runTracking(app, fileKey, preview)
            params = app.Tracker.Parameters(char(fileKey));
            try
                videoReader = VideoReader(char(fileKey));
            catch
                app.appendStatus("Could not read " + fileKey);
                return;
            end

            lastFrame = videoReader.NumFrames;
            if preview
                frameStep = 10;
            else
                frameStep = 1;
            end
            se = strel('disk', app.Tracker.MinRadius);
            [Xgrid, Ygrid] = meshgrid(1:videoReader.Width, 1:videoReader.Height);

            a = params.EyeMajor / 2;
            b = params.EyeMinor / 2;
            theta = -deg2rad(params.EyeOrientation);
            eyeMask = (((Xgrid - params.EyeX) * cos(theta) + (Ygrid - params.EyeY) * sin(theta)).^2 / a^2 + ...
                ((Xgrid - params.EyeX) * sin(theta) - (Ygrid - params.EyeY) * cos(theta)).^2 / b^2) <= 1;
            phiEye = linspace(0, 2 * pi, 50);
            rotationEye = [cos(-theta) sin(-theta); -sin(-theta) cos(-theta)];
            xyEye = rotationEye * [a * cos(phiEye); b * sin(phiEye)];

            nFrames = numel(1:frameStep:lastFrame);
            Frames = NaN(nFrames, 1);
            X = NaN(nFrames, 1);
            Y = NaN(nFrames, 1);
            PupilArea = NaN(nFrames, 1);
            ThresholdUsed = NaN(nFrames, 1);
            MajorAxis = NaN(nFrames, 1);
            MinorAxis = NaN(nFrames, 1);
            Eccentricity = NaN(nFrames, 1);
            Orientation = NaN(nFrames, 1);
            FitQuality = NaN(nFrames, 1);
            EllipseIR = NaN(nFrames, 1);
            BlobIR = NaN(nFrames, 1);
            BoundingBox = NaN(nFrames, 4);

            outIdx = 0;
            app.appendStatus("Tracking " + fileKey);
            for frameNum = 1:frameStep:lastFrame
                outIdx = outIdx + 1;
                Frames(outIdx) = frameNum;
                try
                    frame = read(videoReader, frameNum);
                catch
                    continue;
                end
                if size(frame, 3) == 3
                    frame = rgb2gray(frame);
                end

                binaryImage = (frame > params.Threshold) & eyeMask;
                cleanImage = bwareaopen(imopen(imfill(binaryImage, 'holes'), se), round(params.Area * app.Tracker.MinAreaFrac));
                eyePixels = double(frame(eyeMask));
                meanEye = mean(eyePixels);

                statsEllipse = regionprops(cleanImage, frame, 'Area', 'Centroid', 'MajorAxisLength', 'MinorAxisLength', 'Eccentricity', 'Orientation', 'MeanIntensity');
                if ~isempty(statsEllipse)
                    bestIdx = app.bestRegion(statsEllipse, params.StartX, params.StartY, params.Area);
                    bestEllipse = statsEllipse(bestIdx);
                    X(outIdx) = bestEllipse.Centroid(1);
                    Y(outIdx) = bestEllipse.Centroid(2);
                    PupilArea(outIdx) = bestEllipse.Area;
                    MajorAxis(outIdx) = bestEllipse.MajorAxisLength;
                    MinorAxis(outIdx) = bestEllipse.MinorAxisLength;
                    Eccentricity(outIdx) = bestEllipse.Eccentricity;
                    Orientation(outIdx) = bestEllipse.Orientation;
                    ellipseArea = pi * (bestEllipse.MajorAxisLength / 2) * (bestEllipse.MinorAxisLength / 2);
                    if ellipseArea > 0
                        FitQuality(outIdx) = bestEllipse.Area / ellipseArea;
                    end
                    if meanEye > 0
                        EllipseIR(outIdx) = bestEllipse.MeanIntensity / meanEye;
                    end
                end

                statsBlob = regionprops(cleanImage, frame, 'Area', 'Centroid', 'BoundingBox', 'MeanIntensity');
                if ~isempty(statsBlob)
                    bestBlobIdx = app.bestRegion(statsBlob, params.StartX, params.StartY, params.Area);
                    bestBlob = statsBlob(bestBlobIdx);
                    BoundingBox(outIdx, :) = bestBlob.BoundingBox;
                    if meanEye > 0
                        BlobIR(outIdx) = bestBlob.MeanIntensity / meanEye;
                    end
                end

                if preview
                    row = app.resultRowFromArrays(outIdx, Frames, X, Y, PupilArea, MajorAxis, MinorAxis, Eccentricity, Orientation, BoundingBox, ThresholdUsed, FitQuality, EllipseIR, BlobIR);
                    app.renderPreviewFrame(frame, params, xyEye, row, frameNum, lastFrame);
                    app.IgnorePreviewSlider = true;
                    app.PreviewSlider.Enable = 'on';
                    app.PreviewSlider.Limits = [1 max(2, lastFrame)];
                    app.PreviewSlider.Value = frameNum;
                    app.IgnorePreviewSlider = false;
                    app.PreviewLabel.Text = sprintf('%d / %d', frameNum, lastFrame);
                    drawnow limitrate;
                end
            end

            idx = 1:outIdx;
            resultTable = table(Frames(idx), X(idx), Y(idx), PupilArea(idx), ...
                MajorAxis(idx), MinorAxis(idx), Eccentricity(idx), Orientation(idx), ...
                BoundingBox(idx, :), ThresholdUsed(idx), FitQuality(idx), EllipseIR(idx), BlobIR(idx));
            resultTable.Properties.VariableNames = {'Frame', 'X', 'Y', 'Area', 'MajorAxis', 'MinorAxis', 'Eccentricity', 'Orientation', 'BBox', 'Threshold', 'FitQuality', 'EllipseIR', 'BlobIR'};

            app.Tracker.Results(char(fileKey)) = resultTable;
            outputFile = strrep(char(fileKey), '.mj2', '_pupil.tsv');
            writetable(resultTable, outputFile, 'FileType', 'text', 'Delimiter', '\t');
            app.prepareReviewData(fileKey, resultTable);
            app.appendStatus("Saved tracking results to " + string(outputFile));
        end

        function row = resultRowFromArrays(~, idx, Frames, X, Y, PupilArea, MajorAxis, MinorAxis, Eccentricity, Orientation, BoundingBox, ThresholdUsed, FitQuality, EllipseIR, BlobIR)
            row = table(Frames(idx), X(idx), Y(idx), PupilArea(idx), MajorAxis(idx), MinorAxis(idx), Eccentricity(idx), Orientation(idx), BoundingBox(idx, :), ThresholdUsed(idx), FitQuality(idx), EllipseIR(idx), BlobIR(idx));
            row.Properties.VariableNames = {'Frame', 'X', 'Y', 'Area', 'MajorAxis', 'MinorAxis', 'Eccentricity', 'Orientation', 'BBox', 'Threshold', 'FitQuality', 'EllipseIR', 'BlobIR'};
        end

        function prepareReviewData(app, fileKey, resultTable)
            app.ReviewVideoFile = fileKey;
            app.ReviewResultTable = resultTable;
            app.ReviewParameters = app.Tracker.Parameters(char(fileKey));
            try
                videoReader = VideoReader(char(fileKey));
                app.ReviewLastFrame = videoReader.NumFrames;
            catch
                app.ReviewLastFrame = max(resultTable.Frame);
            end
            a = app.ReviewParameters.EyeMajor / 2;
            b = app.ReviewParameters.EyeMinor / 2;
            theta = -deg2rad(app.ReviewParameters.EyeOrientation);
            phiEye = linspace(0, 2 * pi, 50);
            rotationEye = [cos(-theta) sin(-theta); -sin(-theta) cos(-theta)];
            app.ReviewEyeOutline = rotationEye * [a * cos(phiEye); b * sin(phiEye)];
            app.IgnorePreviewSlider = true;
            app.PreviewSlider.Enable = 'on';
            app.PreviewSlider.Limits = [1 max(2, app.ReviewLastFrame)];
            app.PreviewSlider.Value = min(max(1, resultTable.Frame(end)), app.ReviewLastFrame);
            app.IgnorePreviewSlider = false;
            app.PreviewLabel.Text = sprintf('%d / %d', round(app.PreviewSlider.Value), app.ReviewLastFrame);
            app.WorkspaceTabs.SelectedTab = app.PreviewTab;
            app.onPreviewSliderChanged(app.PreviewSlider.Value);
        end

        function onPreviewSliderChanged(app, value)
            if app.IgnorePreviewSlider || strlength(app.ReviewVideoFile) == 0 || isempty(app.ReviewResultTable)
                return;
            end
            frameNum = max(1, min(app.ReviewLastFrame, round(value)));
            app.IgnorePreviewSlider = true;
            app.PreviewSlider.Value = frameNum;
            app.IgnorePreviewSlider = false;
            try
                videoReader = VideoReader(char(app.ReviewVideoFile));
                frame = read(videoReader, frameNum);
            catch
                return;
            end
            if size(frame, 3) == 3
                frame = rgb2gray(frame);
            end
            [~, rowIdx] = min(abs(app.ReviewResultTable.Frame - frameNum));
            row = app.ReviewResultTable(rowIdx, :);
            app.renderPreviewFrame(frame, app.ReviewParameters, app.ReviewEyeOutline, row, frameNum, app.ReviewLastFrame);
            app.PreviewLabel.Text = sprintf('%d / %d', frameNum, app.ReviewLastFrame);
        end

        function renderPreviewFrame(app, frame, params, xyEye, row, frameNum, lastFrame)
            hold(app.PreviewAxes, 'off');
            imshow(frame, 'Parent', app.PreviewAxes, 'InitialMagnification', 'fit');
            hold(app.PreviewAxes, 'on');
            plot(app.PreviewAxes, xyEye(1, :) + params.EyeX, xyEye(2, :) + params.EyeY, 'b--', 'LineWidth', 1);

            xVal = row.X(1);
            yVal = row.Y(1);
            areaVal = row.Area(1);
            majorAxis = row.MajorAxis(1);
            minorAxis = row.MinorAxis(1);
            orientation = row.Orientation(1);
            fitQuality = row.FitQuality(1);
            ellipseIR = row.EllipseIR(1);
            blobIR = row.BlobIR(1);
            bbox = row.BBox;
            if iscell(bbox)
                bbox = bbox{1};
            end

            if ~isnan(xVal)
                plot(app.PreviewAxes, xVal, yVal, 'g+', 'MarkerSize', 8, 'LineWidth', 2);
                if ~isnan(majorAxis)
                    phiP = linspace(0, 2 * pi, 50);
                    ap = majorAxis / 2;
                    bp = minorAxis / 2;
                    thetaP = pi * orientation / 180;
                    rotationP = [cos(thetaP) sin(thetaP); -sin(thetaP) cos(thetaP)];
                    xyP = rotationP * [ap * cos(phiP); bp * sin(phiP)];
                    plot(app.PreviewAxes, xyP(1, :) + xVal, xyP(2, :) + yVal, 'r-', 'LineWidth', 2);
                end
            end

            if ~isempty(bbox) && ~all(isnan(bbox))
                rectangle(app.PreviewAxes, 'Position', bbox, 'EdgeColor', 'g', 'LineWidth', 2);
            end
            hold(app.PreviewAxes, 'off');

            if isnan(xVal)
                title(app.PreviewAxes, sprintf('Frame %d | BLINK', frameNum), 'Color', 'r');
            else
                title(app.PreviewAxes, sprintf('Frame %d/%d | Area: %d | Fit: %.2f | EIR: %.2f | BIR: %.2f', ...
                    frameNum, lastFrame, areaVal, fitQuality, ellipseIR, blobIR), 'Interpreter', 'none');
            end
        end

        function updateSummary(app)
            if strlength(app.CurrentFile) == 0
                app.SummaryTable.Data = table();
                return;
            end
            rows = table(string(app.CurrentFile), false, false, NaN, NaN, NaN, NaN, ...
                'VariableNames', {'File', 'HasParameters', 'HasResults', 'Threshold', 'StartX', 'StartY', 'FramesTracked'});
            current = char(app.CurrentFile);
            if isKey(app.Tracker.Parameters, current)
                params = app.Tracker.Parameters(current);
                rows.HasParameters(1) = params.Threshold ~= 0;
                rows.Threshold(1) = params.Threshold;
                rows.StartX(1) = params.StartX;
                rows.StartY(1) = params.StartY;
            end
            if isKey(app.Tracker.Results, current)
                resultTable = app.Tracker.Results(current);
                rows.HasResults(1) = ~isempty(resultTable);
                if ~isempty(resultTable)
                    rows.FramesTracked(1) = height(resultTable);
                end
            end
            app.SummaryTable.Data = rows;
        end

        function appendStatus(app, message)
            if isstring(message)
                message = char(message);
            end
            if isempty(app.StatusTextArea.Value)
                app.StatusTextArea.Value = {message};
            else
                app.StatusTextArea.Value = [app.StatusTextArea.Value; {message}];
            end
            drawnow limitrate;
        end

        function ok = ensureVideoLoaded(app)
            ok = false;
            if strlength(app.CurrentFile) == 0
                uialert(app.UIFigure, 'Select a file first.', 'No File Selected');
                return;
            end
            if strlength(app.CurrentVideoFile) > 0 && app.CurrentVideoFile == app.CurrentFile && ~isempty(app.CurrentVideoReader)
                ok = true;
                return;
            end
            try
                app.CurrentVideoReader = VideoReader(char(app.CurrentFile));
                app.CurrentVideoFile = app.CurrentFile;
                ok = true;
            catch exc
                uialert(app.UIFigure, exc.message, 'Video Read Error');
            end
        end

        function framePool = getFramePool(~, videoReader, n)
            totalFrames = round(videoReader.Duration * videoReader.FrameRate);
            n = min(n, totalFrames);
            frameIdx = sort(randperm(totalFrames, n));
            framePool = zeros(videoReader.Height, videoReader.Width, n, sprintf('uint%d', videoReader.BitsPerPixel));
            for idx = 1:n
                try
                    framePool(:, :, idx) = read(videoReader, frameIdx(idx));
                catch
                end
            end
        end

        function idx = bestRegion(~, stats, refX, refY, targetArea)
            n = numel(stats);
            if n == 1
                idx = 1;
                return;
            end
            centroids = vertcat(stats.Centroid);
            dists = sqrt((centroids(:, 1) - refX).^2 + (centroids(:, 2) - refY).^2);
            areas = [stats.Area]';
            intensities = [stats.MeanIntensity]';
            [~, order] = sort(dists, 'ascend');
            rankPos = zeros(n, 1);
            rankPos(order) = 1:n;
            [~, order] = sort(abs(areas - targetArea), 'ascend');
            rankArea = zeros(n, 1);
            rankArea(order) = 1:n;
            [~, order] = sort(intensities, 'descend');
            rankIntensity = zeros(n, 1);
            rankIntensity(order) = 1:n;
            [~, idx] = min(rankPos + rankArea + rankIntensity);
        end
    end
end
%{
classdef PupilTrackerApp < matlab.apps.AppBase
    % sbx.PupilTrackerApp  App Designer-style interface for pupil setup and tracking.

    properties (Access = public)
        UIFigure
        MainGrid
        SidebarGrid
        WorkspaceTabs
        FileListLabel
        FileListBox
        InitializeButton
        NewAverageButton
        ConfirmEyeButton
        SaveParametersButton
        RestartButton
        SkipButton
        TrackCurrentButton
        TrackCurrentPreviewButton
        TrackAllButton
        TrackAllPreviewButton
        StatusTextArea
        EyeTab
        PupilTab
        PreviewTab
        SummaryTab
        EyeAxes
        PupilGrid
        PreviewAxes
        PreviewSlider
        PreviewLabel
        SummaryTable
    end

    properties (Access = private)
        Tracker
        FileKeys = strings(0,1)
        CurrentFile = ""
        CurrentVideoReader = []
        CurrentVideoFile = ""
        AverageFrame = []
        EyeRoi = []
        PupilRois = gobjects(0)
        PupilAxes = gobjects(0)
        SampleFramePool = []
        PupilFramePool = []
        QuantPositions = []
        EyeMask = []
        EyeCenter = []
        EyeSemiAxes = []
        EyeAngle = 0
        ReviewVideoFile = ""
        ReviewResultTable = table()
        ReviewParameters = table()
        ReviewLastFrame = 0
        ReviewEyeOutline = zeros(2, 0)
        IgnorePreviewSlider = false
    end

    methods (Access = public)
        function app = PupilTrackerApp(tracker, mode)
            arguments
                tracker (1,1) sbx.PupilTracker
                mode (1,1) string {mustBeMember(mode, ["initialize", "track"])} = "initialize"
            end

            app.Tracker = tracker;
            app.createComponents();
            registerApp(app, app.UIFigure);
            app.populateFileList();
            app.switchMode(mode);

            if ~isempty(app.FileKeys)
                app.selectFile(app.FileKeys(1));
            else
                app.appendStatus('No files are currently queued in this PupilTracker object.');
            end

            if nargout == 0
                clear app
            end
        end

        function delete(app)
            if ~isempty(app.UIFigure) && isvalid(app.UIFigure)
                delete(app.UIFigure);
            end
        end

        function switchMode(app, mode)
            if mode == "track"
                app.WorkspaceTabs.SelectedTab = app.PreviewTab;
                app.appendStatus('Track mode selected. Use Track Current or Track All.');
            else
                app.WorkspaceTabs.SelectedTab = app.EyeTab;
                app.appendStatus('Initialize mode selected. Load an average frame, draw the eye ROI, then save pupil parameters.');
            end
        end
    end

    methods (Access = private)
        function createComponents(app)
            app.UIFigure = uifigure('Name', 'Pupil Tracker App', 'Position', [100 100 1500 920]);
            app.UIFigure.CloseRequestFcn = @(~,~) delete(app);

            app.MainGrid = uigridlayout(app.UIFigure, [1 2]);
            app.MainGrid.ColumnWidth = {280, '1x'};
            app.MainGrid.RowHeight = {'1x'};
            app.MainGrid.Padding = [8 8 8 8];

            app.SidebarGrid = uigridlayout(app.MainGrid, [14 1]);
            app.SidebarGrid.Layout.Row = 1;
            app.SidebarGrid.Layout.Column = 1;
            app.SidebarGrid.RowHeight = {22, 170, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, '1x', 180};
            app.SidebarGrid.ColumnWidth = {'1x'};

            app.FileListLabel = uilabel(app.SidebarGrid, 'Text', 'Queued Files');
            app.FileListLabel.Layout.Row = 1;

            app.FileListBox = uilistbox(app.SidebarGrid, 'ValueChangedFcn', @(~,~) app.onFileChanged());
            app.FileListBox.Layout.Row = 2;

            app.InitializeButton = uibutton(app.SidebarGrid, 'push', 'Text', 'Load Eye Average', ...
                'ButtonPushedFcn', @(~,~) app.loadEyeAverage());
            app.InitializeButton.Layout.Row = 3;

            app.NewAverageButton = uibutton(app.SidebarGrid, 'push', 'Text', 'New Average', ...
                'ButtonPushedFcn', @(~,~) app.loadEyeAverage());
            app.NewAverageButton.Layout.Row = 4;

            app.ConfirmEyeButton = uibutton(app.SidebarGrid, 'push', 'Text', 'Confirm Eye / Load Pupils', ...
                'ButtonPushedFcn', @(~,~) app.confirmEye());
            app.ConfirmEyeButton.Layout.Row = 5;

            app.SaveParametersButton = uibutton(app.SidebarGrid, 'push', 'Text', 'Save Parameters', ...
                'ButtonPushedFcn', @(~,~) app.saveParameters());
            app.SaveParametersButton.Layout.Row = 6;

            app.RestartButton = uibutton(app.SidebarGrid, 'push', 'Text', 'Restart Current File', ...
                'ButtonPushedFcn', @(~,~) app.restartCurrentFile());
            app.RestartButton.Layout.Row = 7;

            app.SkipButton = uibutton(app.SidebarGrid, 'push', 'Text', 'Skip Current File', ...
                'ButtonPushedFcn', @(~,~) app.skipCurrentFile());
            app.SkipButton.Layout.Row = 8;

            app.TrackCurrentPreviewButton = uibutton(app.SidebarGrid, 'push', 'Text', 'Track Current Preview', ...
                'ButtonPushedFcn', @(~,~) app.trackFiles(false, true));
            app.TrackCurrentPreviewButton.Layout.Row = 9;

            app.TrackCurrentButton = uibutton(app.SidebarGrid, 'push', 'Text', 'Track Current', ...
                'ButtonPushedFcn', @(~,~) app.trackFiles(false, false));
            app.TrackCurrentButton.Layout.Row = 10;

            app.TrackAllPreviewButton = uibutton(app.SidebarGrid, 'push', 'Text', 'Track All Preview', ...
                'ButtonPushedFcn', @(~,~) app.trackFiles(true, true));
            app.TrackAllPreviewButton.Layout.Row = 11;

            app.TrackAllButton = uibutton(app.SidebarGrid, 'push', 'Text', 'Track All', ...
                'ButtonPushedFcn', @(~,~) app.trackFiles(true, false));
            app.TrackAllButton.Layout.Row = 12;

            app.StatusTextArea = uitextarea(app.SidebarGrid, 'Editable', 'off');
            app.StatusTextArea.Layout.Row = 14;
            app.StatusTextArea.Value = {'Ready.'};

            app.WorkspaceTabs = uitabgroup(app.MainGrid);
            app.WorkspaceTabs.Layout.Row = 1;
            app.WorkspaceTabs.Layout.Column = 2;

            app.EyeTab = uitab(app.WorkspaceTabs, 'Title', 'Eye ROI');
            app.EyeAxes = uiaxes(app.EyeTab, 'Position', [10 10 1180 820]);
            app.configureImageAxes(app.EyeAxes);

            app.PupilTab = uitab(app.WorkspaceTabs, 'Title', 'Pupil ROIs');
            app.PupilGrid = uigridlayout(app.PupilTab, [2 3]);
            app.PupilGrid.RowHeight = {'1x', '1x'};
            app.PupilGrid.ColumnWidth = {'1x', '1x', '1x'};
            app.PupilGrid.Padding = [8 8 8 8];
            app.PupilGrid.RowSpacing = 8;
            app.PupilGrid.ColumnSpacing = 8;
            app.PupilAxes = gobjects(app.Tracker.NrPupilImages, 1);
            for idx = 1:app.Tracker.NrPupilImages
                app.PupilAxes(idx) = uiaxes(app.PupilGrid);
                app.configureImageAxes(app.PupilAxes(idx));
            end

            app.PreviewTab = uitab(app.WorkspaceTabs, 'Title', 'Tracking Preview');
            previewGrid = uigridlayout(app.PreviewTab, [3 1]);
            previewGrid.RowHeight = {'1x', 40, 30};
            previewGrid.ColumnWidth = {'1x'};
            previewGrid.Padding = [8 8 8 8];
            previewGrid.RowSpacing = 8;

            app.PreviewAxes = uiaxes(previewGrid);
            app.PreviewAxes.Layout.Row = 1;
            app.configureImageAxes(app.PreviewAxes);

            app.PreviewSlider = uislider(previewGrid, ...
                'Limits', [1 2], 'Value', 1, 'Enable', 'off', ...
                'ValueChangedFcn', @(src,~) app.onPreviewSliderChanged(src.Value), ...
                'ValueChangingFcn', @(~,evt) app.onPreviewSliderChanged(evt.Value));
            app.PreviewSlider.Layout.Row = 2;

            app.PreviewLabel = uilabel(previewGrid, 'Text', '0 / 0');
            app.PreviewLabel.Layout.Row = 3;

            app.SummaryTab = uitab(app.WorkspaceTabs, 'Title', 'Summary');
            app.SummaryTable = uitable(app.SummaryTab, 'Position', [10 10 1180 820]);
            app.SummaryTable.ColumnEditable = false;
        end

        function configureImageAxes(~, ax)
            ax.Toolbar.Visible = 'off';
            ax.XTick = [];
            ax.YTick = [];
            ax.Box = 'on';
            ax.YDir = 'reverse';
            ax.DataAspectRatio = [1 1 1];
        end

        function populateFileList(app)
            app.FileKeys = string(keys(app.Tracker.Parameters))';
            if isempty(app.FileKeys)
                app.FileListBox.Items = {};
                app.FileListBox.ItemsData = {};
                return;
            end

            labels = cell(numel(app.FileKeys), 1);
            itemsData = cell(numel(app.FileKeys), 1);
            for idx = 1:numel(app.FileKeys)
                fileKey = app.FileKeys(idx);
                labels{idx} = app.makeFileLabel(fileKey);
                itemsData{idx} = char(fileKey);
            end
            app.FileListBox.Items = labels;
            app.FileListBox.ItemsData = itemsData;
            if isempty(app.FileListBox.Value) || ~any(strcmp(app.FileListBox.Value, itemsData))
                app.FileListBox.Value = itemsData{1};
            end
        end

        function label = makeFileLabel(app, fileKey)
            [~, name, ext] = fileparts(char(fileKey));
            hasParams = isKey(app.Tracker.Parameters, char(fileKey)) && app.Tracker.Parameters(char(fileKey)).Threshold ~= 0;
            hasResults = isKey(app.Tracker.Results, char(fileKey)) && ~isempty(app.Tracker.Results(char(fileKey)));
            if hasParams
                pFlag = 'P';
            else
                pFlag = '-';
            end
            if hasResults
                rFlag = 'R';
            else
                rFlag = '-';
            end
            label = sprintf('[%s%s] %s%s', pFlag, rFlag, name, ext);
        end

        function onFileChanged(app)
            if isempty(app.FileListBox.ItemsData)
                return;
            end
            app.selectFile(string(app.FileListBox.Value));
        end

        function selectFile(app, fileKey)
            if strlength(fileKey) == 0
                return;
            end
            app.CurrentFile = fileKey;
            app.resetTransientState();
            app.updateSummary();
            app.appendStatus("Selected " + fileKey);

            if isKey(app.Tracker.Results, char(fileKey)) && ~isempty(app.Tracker.Results(char(fileKey)))
                app.prepareReviewData(fileKey, app.Tracker.Results(char(fileKey)));
            else
                app.PreviewSlider.Enable = 'off';
                app.PreviewLabel.Text = '0 / 0';
                hold(app.PreviewAxes, 'off');
                cla(app.PreviewAxes);
            end
        end

        function resetTransientState(app)
            if ~isempty(app.EyeRoi) && isvalid(app.EyeRoi)
                delete(app.EyeRoi);
            end
            app.EyeRoi = [];
            if ~isempty(app.PupilRois)
                for idx = 1:numel(app.PupilRois)
                    if isgraphics(app.PupilRois(idx))
                        delete(app.PupilRois(idx));
                    end
                end
            end
            app.PupilRois = gobjects(0);
            app.SampleFramePool = [];
            app.PupilFramePool = [];
            app.QuantPositions = [];
            app.EyeMask = [];
            app.EyeCenter = [];
            app.EyeSemiAxes = [];
            app.EyeAngle = 0;
            hold(app.EyeAxes, 'off');
            cla(app.EyeAxes);
            for idx = 1:numel(app.PupilAxes)
                hold(app.PupilAxes(idx), 'off');
                cla(app.PupilAxes(idx));
            end
        end

        function loadEyeAverage(app)
            if ~app.ensureVideoLoaded()
                return;
            end
            app.SampleFramePool = app.getFramePool(app.CurrentVideoReader, app.Tracker.MaxSampleFrames);
            brightnessPerFrame = squeeze(mean(app.SampleFramePool, [1 2]));
            keep = brightnessPerFrame > prctile(brightnessPerFrame, 50);
            if ~any(keep)
                keep = true(size(brightnessPerFrame));
            end
            app.AverageFrame = uint8(mean(double(app.SampleFramePool(:, :, keep)), 3));

            hold(app.EyeAxes, 'off');
            imshow(app.AverageFrame, 'Parent', app.EyeAxes, 'InitialMagnification', 'fit');
            title(app.EyeAxes, 'Draw eye boundary, then click Confirm Eye / Load Pupils.');
            if ~isempty(app.EyeRoi) && isvalid(app.EyeRoi)
                delete(app.EyeRoi);
            end
            app.EyeRoi = drawellipse(app.EyeAxes, 'Color', 'b', 'FaceAlpha', 0.15);
            app.WorkspaceTabs.SelectedTab = app.EyeTab;
            app.appendStatus('Loaded a new average eye frame.');
        end

        function confirmEye(app)
            if isempty(app.EyeRoi) || ~isvalid(app.EyeRoi)
                uialert(app.UIFigure, 'Draw an eye ellipse first.', 'Eye ROI Required');
                return;
            end
            if isempty(app.SampleFramePool)
                uialert(app.UIFigure, 'Load an average eye frame first.', 'Sample Frames Required');
                return;
            end

            app.EyeCenter = app.EyeRoi.Center;
            app.EyeSemiAxes = app.EyeRoi.SemiAxes;
            app.EyeAngle = app.EyeRoi.RotationAngle;

            [Xg, Yg] = meshgrid(1:app.CurrentVideoReader.Width, 1:app.CurrentVideoReader.Height);
            aEye = app.EyeSemiAxes(1);
            bEye = app.EyeSemiAxes(2);
            thetaEye = -deg2rad(app.EyeAngle);
            app.EyeMask = (((Xg - app.EyeCenter(1)) * cos(thetaEye) + (Yg - app.EyeCenter(2)) * sin(thetaEye)).^2 / aEye^2 + ...
                ((Xg - app.EyeCenter(1)) * sin(thetaEye) - (Yg - app.EyeCenter(2)) * cos(thetaEye)).^2 / bEye^2) <= 1;

            eyePixelsRef = double(app.AverageFrame(app.EyeMask));
            eyeFloorRef = prctile(eyePixelsRef, 5);
            nPool = size(app.SampleFramePool, 3);
            eyeBrightness = zeros(nPool, 1);
            keepFrame = true(nPool, 1);
            for frameIdx = 1:nPool
                frameSlice = double(app.SampleFramePool(:, :, frameIdx));
                eyeBrightness(frameIdx) = mean(frameSlice(app.EyeMask));
                if eyeBrightness(frameIdx) < eyeFloorRef
                    keepFrame(frameIdx) = false;
                end
            end
            keptIdx = find(keepFrame);
            if numel(keptIdx) < app.Tracker.NrPupilImages
                keptIdx = 1:nPool;
            end
            [~, rankAmongKept] = sort(eyeBrightness(keptIdx));
            rankOrder = keptIdx(rankAmongKept);
            nKept = numel(rankOrder);
            app.QuantPositions = round(linspace(1, nKept, app.Tracker.NrPupilImages));
            tileIdx = rankOrder(app.QuantPositions);
            app.PupilFramePool = app.SampleFramePool(:, :, tileIdx);

            phiEye = linspace(0, 2 * pi, 100);
            alpha = deg2rad(app.EyeAngle);
            rotationEye = [cos(alpha) sin(alpha); -sin(alpha) cos(alpha)];
            xyEye = rotationEye * [aEye * cos(phiEye); bEye * sin(phiEye)];

            if ~isempty(app.PupilRois)
                for idx = 1:numel(app.PupilRois)
                    if isgraphics(app.PupilRois(idx))
                        delete(app.PupilRois(idx));
                    end
                end
            end
            app.PupilRois = gobjects(app.Tracker.NrPupilImages, 1);
            for idx = 1:app.Tracker.NrPupilImages
                ax = app.PupilAxes(idx);
                hold(ax, 'off');
                imshow(app.PupilFramePool(:, :, idx), 'Parent', ax, 'InitialMagnification', 'fit');
                hold(ax, 'on');
                plot(ax, xyEye(1, :) + app.EyeCenter(1), xyEye(2, :) + app.EyeCenter(2), 'b-', 'LineWidth', 1.5);
                title(ax, sprintf('Rank %d/%d', app.QuantPositions(idx), nKept));
                app.PupilRois(idx) = drawellipse(ax, 'Color', 'r', 'FaceAlpha', 0.15);
                hold(ax, 'off');
            end
            app.WorkspaceTabs.SelectedTab = app.PupilTab;
            app.appendStatus('Eye ROI confirmed. Draw or adjust the red pupil ellipses, then click Save Parameters.');
        end

        function saveParameters(app)
            if strlength(app.CurrentFile) == 0 || isempty(app.PupilFramePool) || isempty(app.EyeMask)
                uialert(app.UIFigure, 'Complete eye and pupil selection first.', 'Parameters Not Ready');
                return;
            end
            allPupilPixels = [];
            allEyePixels = [];
            validCentroids = zeros(0, 2);
            validAreas = zeros(0, 1);
            for idx = 1:numel(app.PupilRois)
                if ~isgraphics(app.PupilRois(idx))
                    continue;
                end
                cx = app.PupilRois(idx).Center(1);
                cy = app.PupilRois(idx).Center(2);
                row = round(cy);
                col = round(cx);
                if row < 1 || row > app.CurrentVideoReader.Height || col < 1 || col > app.CurrentVideoReader.Width
                    continue;
                end
                if ~app.EyeMask(row, col)
                    continue;
                end
                mask = createMask(app.PupilRois(idx));
                frame = app.PupilFramePool(:, :, idx);
                thisEyePixels = double(frame(app.EyeMask & ~mask));
                thisPupilPixels = double(frame(mask));
                allPupilPixels = [allPupilPixels; thisPupilPixels]; %#ok<AGROW>
                allEyePixels = [allEyePixels; thisEyePixels]; %#ok<AGROW>
                validCentroids(end+1, :) = [cx, cy]; %#ok<AGROW>
                validAreas(end+1) = pi * prod(app.PupilRois(idx).SemiAxes); %#ok<AGROW>
            end

            if isempty(allPupilPixels)
                uialert(app.UIFigure, 'No valid pupil ellipses are inside the eye ROI.', 'No Valid Pupil ROIs');
                return;
            end

            app.Tracker.Parameters(char(app.CurrentFile)) = table( ...
                (median(allPupilPixels) + median(allEyePixels)) / 2, ...
                median(allEyePixels), ...
                median(validCentroids(:, 1)), median(validCentroids(:, 2)), ...
                median(validAreas), ...
                app.EyeCenter(1), app.EyeCenter(2), ...
                app.EyeSemiAxes(1) * 2, app.EyeSemiAxes(2) * 2, app.EyeAngle, ...
                'VariableNames', {'Threshold', 'Floor', 'StartX', 'StartY', 'Area', 'EyeX', 'EyeY', 'EyeMajor', 'EyeMinor', 'EyeOrientation'});

            [folderName, fileName, ext] = fileparts(char(app.CurrentFile));
            paramStruct = table2struct(app.Tracker.Parameters(char(app.CurrentFile)));
            paramStruct.FileName = [fileName ext];
            jsonText = jsonencode(paramStruct, 'PrettyPrint', true);
            jsonFile = fullfile(folderName, strrep([fileName ext], '.mj2', '_pupil.json'));
            fileId = fopen(jsonFile, 'w');
            if fileId == -1
                uialert(app.UIFigure, ['Could not write ' jsonFile], 'Write Error');
                return;
            end
            fwrite(fileId, jsonText, 'char');
            fclose(fileId);

            app.populateFileList();
            app.updateSummary();
            app.appendStatus("Saved parameters to " + jsonFile);
        end

        function restartCurrentFile(app)
            if strlength(app.CurrentFile) == 0
                return;
            end
            app.resetTransientState();
            app.loadEyeAverage();
            app.appendStatus('Restarted current file initialization.');
        end

        function skipCurrentFile(app)
            if strlength(app.CurrentFile) == 0
                return;
            end
            current = char(app.CurrentFile);
            if isKey(app.Tracker.Parameters, current)
                remove(app.Tracker.Parameters, current);
            end
            if isKey(app.Tracker.Results, current)
                remove(app.Tracker.Results, current);
            end
            app.appendStatus("Skipped " + app.CurrentFile);
            app.populateFileList();
            if ~isempty(app.FileKeys)
                app.selectFile(app.FileKeys(1));
            else
                app.CurrentFile = "";
                app.resetTransientState();
                app.SummaryTable.Data = table();
            end
        end

        function trackFiles(app, allFiles, preview)
            if isempty(app.FileKeys)
                uialert(app.UIFigure, 'There are no queued files to track.', 'Nothing To Track');
                return;
            end
            if allFiles
                fileSet = app.FileKeys;
            else
                fileSet = app.CurrentFile;
            end
            fileSet = fileSet(strlength(fileSet) > 0);
            if isempty(fileSet)
                return;
            end
            app.WorkspaceTabs.SelectedTab = app.PreviewTab;
            for idx = 1:numel(fileSet)
                fileKey = fileSet(idx);
                if ~isKey(app.Tracker.Parameters, char(fileKey)) || app.Tracker.Parameters(char(fileKey)).Threshold == 0
                    app.appendStatus("Skipping " + fileKey + ' because parameters are not initialized.');
                    continue;
                end
                app.runTracking(fileKey, preview);
                app.populateFileList();
            end
            app.updateSummary();
        end

        function runTracking(app, fileKey, preview)
            params = app.Tracker.Parameters(char(fileKey));
            try
                videoReader = VideoReader(char(fileKey));
            catch
                app.appendStatus("Could not read " + fileKey);
                return;
            end

            lastFrame = videoReader.NumFrames;
            if preview
                frameStep = 10;
            else
                frameStep = 1;
            end
            se = strel('disk', app.Tracker.MinRadius);
            [Xgrid, Ygrid] = meshgrid(1:videoReader.Width, 1:videoReader.Height);

            a = params.EyeMajor / 2;
            b = params.EyeMinor / 2;
            theta = -deg2rad(params.EyeOrientation);
            eyeMask = (((Xgrid - params.EyeX) * cos(theta) + (Ygrid - params.EyeY) * sin(theta)).^2 / a^2 + ...
                ((Xgrid - params.EyeX) * sin(theta) - (Ygrid - params.EyeY) * cos(theta)).^2 / b^2) <= 1;
            phiEye = linspace(0, 2 * pi, 50);
            rotationEye = [cos(-theta) sin(-theta); -sin(-theta) cos(-theta)];
            xyEye = rotationEye * [a * cos(phiEye); b * sin(phiEye)];

            nFrames = numel(1:frameStep:lastFrame);
            Frames = NaN(nFrames, 1);
            X = NaN(nFrames, 1);
            Y = NaN(nFrames, 1);
            PupilArea = NaN(nFrames, 1);
            ThresholdUsed = NaN(nFrames, 1);
            MajorAxis = NaN(nFrames, 1);
            MinorAxis = NaN(nFrames, 1);
            Eccentricity = NaN(nFrames, 1);
            Orientation = NaN(nFrames, 1);
            FitQuality = NaN(nFrames, 1);
            EllipseIR = NaN(nFrames, 1);
            BlobIR = NaN(nFrames, 1);
            BoundingBox = NaN(nFrames, 4);

            outIdx = 0;
            app.appendStatus("Tracking " + fileKey);
            for frameNum = 1:frameStep:lastFrame
                outIdx = outIdx + 1;
                Frames(outIdx) = frameNum;
                try
                    frame = read(videoReader, frameNum);
                catch
                    continue;
                end
                if size(frame, 3) == 3
                    frame = rgb2gray(frame);
                end

                binaryImage = (frame > params.Threshold) & eyeMask;
                cleanImage = bwareaopen(imopen(imfill(binaryImage, 'holes'), se), round(params.Area * app.Tracker.MinAreaFrac));
                eyePixels = double(frame(eyeMask));
                meanEye = mean(eyePixels);

                statsEllipse = regionprops(cleanImage, frame, 'Area', 'Centroid', 'MajorAxisLength', 'MinorAxisLength', 'Eccentricity', 'Orientation', 'MeanIntensity');
                if ~isempty(statsEllipse)
                    bestIdx = app.bestRegion(statsEllipse, params.StartX, params.StartY, params.Area);
                    bestEllipse = statsEllipse(bestIdx);
                    X(outIdx) = bestEllipse.Centroid(1);
                    Y(outIdx) = bestEllipse.Centroid(2);
                    PupilArea(outIdx) = bestEllipse.Area;
                    MajorAxis(outIdx) = bestEllipse.MajorAxisLength;
                    MinorAxis(outIdx) = bestEllipse.MinorAxisLength;
                    Eccentricity(outIdx) = bestEllipse.Eccentricity;
                    Orientation(outIdx) = bestEllipse.Orientation;
                    ellipseArea = pi * (bestEllipse.MajorAxisLength / 2) * (bestEllipse.MinorAxisLength / 2);
                    if ellipseArea > 0
                        FitQuality(outIdx) = bestEllipse.Area / ellipseArea;
                    end
                    if meanEye > 0
                        EllipseIR(outIdx) = bestEllipse.MeanIntensity / meanEye;
                    end
                end

                statsBlob = regionprops(cleanImage, frame, 'Area', 'Centroid', 'BoundingBox', 'MeanIntensity');
                if ~isempty(statsBlob)
                    bestBlobIdx = app.bestRegion(statsBlob, params.StartX, params.StartY, params.Area);
                    bestBlob = statsBlob(bestBlobIdx);
                    BoundingBox(outIdx, :) = bestBlob.BoundingBox;
                    if meanEye > 0
                        BlobIR(outIdx) = bestBlob.MeanIntensity / meanEye;
                    end
                end

                if preview
                    row = app.resultRowFromArrays(outIdx, Frames, X, Y, PupilArea, MajorAxis, MinorAxis, Eccentricity, Orientation, BoundingBox, ThresholdUsed, FitQuality, EllipseIR, BlobIR);
                    app.renderPreviewFrame(frame, params, xyEye, row, frameNum, lastFrame);
                    app.IgnorePreviewSlider = true;
                    app.PreviewSlider.Enable = 'on';
                    app.PreviewSlider.Limits = [1 max(2, lastFrame)];
                    app.PreviewSlider.Value = frameNum;
                    app.IgnorePreviewSlider = false;
                    app.PreviewLabel.Text = sprintf('%d / %d', frameNum, lastFrame);
                    drawnow limitrate;
                end
            end

            idx = 1:outIdx;
            resultTable = table(Frames(idx), X(idx), Y(idx), PupilArea(idx), ...
                MajorAxis(idx), MinorAxis(idx), Eccentricity(idx), Orientation(idx), ...
                BoundingBox(idx, :), ThresholdUsed(idx), FitQuality(idx), EllipseIR(idx), BlobIR(idx));
            resultTable.Properties.VariableNames = {'Frame', 'X', 'Y', 'Area', 'MajorAxis', 'MinorAxis', 'Eccentricity', 'Orientation', 'BBox', 'Threshold', 'FitQuality', 'EllipseIR', 'BlobIR'};

            app.Tracker.Results(char(fileKey)) = resultTable;
            outputFile = strrep(char(fileKey), '.mj2', '_pupil.tsv');
            writetable(resultTable, outputFile, 'FileType', 'text', 'Delimiter', '\t');
            app.prepareReviewData(fileKey, resultTable);
            app.appendStatus("Saved tracking results to " + string(outputFile));
        end

        function row = resultRowFromArrays(~, idx, Frames, X, Y, PupilArea, MajorAxis, MinorAxis, Eccentricity, Orientation, BoundingBox, ThresholdUsed, FitQuality, EllipseIR, BlobIR)
            row = table(Frames(idx), X(idx), Y(idx), PupilArea(idx), MajorAxis(idx), MinorAxis(idx), Eccentricity(idx), Orientation(idx), BoundingBox(idx, :), ThresholdUsed(idx), FitQuality(idx), EllipseIR(idx), BlobIR(idx));
            row.Properties.VariableNames = {'Frame', 'X', 'Y', 'Area', 'MajorAxis', 'MinorAxis', 'Eccentricity', 'Orientation', 'BBox', 'Threshold', 'FitQuality', 'EllipseIR', 'BlobIR'};
        end

        function prepareReviewData(app, fileKey, resultTable)
            app.ReviewVideoFile = fileKey;
            app.ReviewResultTable = resultTable;
            app.ReviewParameters = app.Tracker.Parameters(char(fileKey));
            try
                videoReader = VideoReader(char(fileKey));
                app.ReviewLastFrame = videoReader.NumFrames;
            catch
                app.ReviewLastFrame = max(resultTable.Frame);
            end
            a = app.ReviewParameters.EyeMajor / 2;
            b = app.ReviewParameters.EyeMinor / 2;
            theta = -deg2rad(app.ReviewParameters.EyeOrientation);
            phiEye = linspace(0, 2 * pi, 50);
            rotationEye = [cos(-theta) sin(-theta); -sin(-theta) cos(-theta)];
            app.ReviewEyeOutline = rotationEye * [a * cos(phiEye); b * sin(phiEye)];
            app.IgnorePreviewSlider = true;
            app.PreviewSlider.Enable = 'on';
            app.PreviewSlider.Limits = [1 max(2, app.ReviewLastFrame)];
            app.PreviewSlider.Value = min(max(1, resultTable.Frame(end)), app.ReviewLastFrame);
            app.IgnorePreviewSlider = false;
            app.PreviewLabel.Text = sprintf('%d / %d', round(app.PreviewSlider.Value), app.ReviewLastFrame);
            app.WorkspaceTabs.SelectedTab = app.PreviewTab;
            app.onPreviewSliderChanged(app.PreviewSlider.Value);
        end

        function onPreviewSliderChanged(app, value)
            if app.IgnorePreviewSlider || strlength(app.ReviewVideoFile) == 0 || isempty(app.ReviewResultTable)
                return;
            end
            frameNum = max(1, min(app.ReviewLastFrame, round(value)));
            app.IgnorePreviewSlider = true;
            app.PreviewSlider.Value = frameNum;
            app.IgnorePreviewSlider = false;
            try
                videoReader = VideoReader(char(app.ReviewVideoFile));
                frame = read(videoReader, frameNum);
            catch
                return;
            end
            if size(frame, 3) == 3
                frame = rgb2gray(frame);
            end
            [~, rowIdx] = min(abs(app.ReviewResultTable.Frame - frameNum));
            row = app.ReviewResultTable(rowIdx, :);
            app.renderPreviewFrame(frame, app.ReviewParameters, app.ReviewEyeOutline, row, frameNum, app.ReviewLastFrame);
            app.PreviewLabel.Text = sprintf('%d / %d', frameNum, app.ReviewLastFrame);
        end

        function renderPreviewFrame(app, frame, params, xyEye, row, frameNum, lastFrame)
            hold(app.PreviewAxes, 'off');
            imshow(frame, 'Parent', app.PreviewAxes, 'InitialMagnification', 'fit');
            hold(app.PreviewAxes, 'on');
            plot(app.PreviewAxes, xyEye(1, :) + params.EyeX, xyEye(2, :) + params.EyeY, 'b--', 'LineWidth', 1);
            xVal = row.X(1);
            yVal = row.Y(1);
            areaVal = row.Area(1);
            majorAxis = row.MajorAxis(1);
            minorAxis = row.MinorAxis(1);
            orientation = row.Orientation(1);
            fitQuality = row.FitQuality(1);
            ellipseIR = row.EllipseIR(1);
            blobIR = row.BlobIR(1);
            bbox = row.BBox;
            if iscell(bbox)
                bbox = bbox{1};
            end

            if ~isnan(xVal)
                plot(app.PreviewAxes, xVal, yVal, 'g+', 'MarkerSize', 8, 'LineWidth', 2);
                if ~isnan(majorAxis)
                    phiP = linspace(0, 2 * pi, 50);
                    ap = majorAxis / 2;
                    bp = minorAxis / 2;
                    thetaP = pi * orientation / 180;
                    rotationP = [cos(thetaP) sin(thetaP); -sin(thetaP) cos(thetaP)];
                    xyP = rotationP * [ap * cos(phiP); bp * sin(phiP)];
                    plot(app.PreviewAxes, xyP(1, :) + xVal, xyP(2, :) + yVal, 'r-', 'LineWidth', 2);
                end
            end

            if ~isempty(bbox) && ~all(isnan(bbox))
                rectangle(app.PreviewAxes, 'Position', bbox, 'EdgeColor', 'g', 'LineWidth', 2);
            end
            hold(app.PreviewAxes, 'off');
            if isnan(xVal)
                title(app.PreviewAxes, sprintf('Frame %d | BLINK', frameNum), 'Color', 'r');
            else
                title(app.PreviewAxes, sprintf('Frame %d/%d | Area: %d | Fit: %.2f | EIR: %.2f | BIR: %.2f', ...
                    frameNum, lastFrame, areaVal, fitQuality, ellipseIR, blobIR), 'Interpreter', 'none');
            end
        end

        function updateSummary(app)
            if strlength(app.CurrentFile) == 0
                app.SummaryTable.Data = table();
                return;
            end
            rows = table(string(app.CurrentFile), false, false, NaN, NaN, NaN, NaN, ...
                'VariableNames', {'File', 'HasParameters', 'HasResults', 'Threshold', 'StartX', 'StartY', 'FramesTracked'});
            current = char(app.CurrentFile);
            if isKey(app.Tracker.Parameters, current)
                params = app.Tracker.Parameters(current);
                rows.HasParameters(1) = params.Threshold ~= 0;
                rows.Threshold(1) = params.Threshold;
                rows.StartX(1) = params.StartX;
                rows.StartY(1) = params.StartY;
            end
            if isKey(app.Tracker.Results, current)
                resultTable = app.Tracker.Results(current);
                rows.HasResults(1) = ~isempty(resultTable);
                if ~isempty(resultTable)
                    rows.FramesTracked(1) = height(resultTable);
                end
            end
            app.SummaryTable.Data = rows;
        end

        function appendStatus(app, message)
            if isstring(message)
                message = char(message);
            end
            if isempty(app.StatusTextArea.Value)
                app.StatusTextArea.Value = {message};
            else
                app.StatusTextArea.Value = [app.StatusTextArea.Value; {message}];
            end
            drawnow limitrate;
        end

        function ok = ensureVideoLoaded(app)
            ok = false;
            if strlength(app.CurrentFile) == 0
                uialert(app.UIFigure, 'Select a file first.', 'No File Selected');
                return;
            end
            if strlength(app.CurrentVideoFile) > 0 && app.CurrentVideoFile == app.CurrentFile && ~isempty(app.CurrentVideoReader)
                ok = true;
                return;
            end
            try
                app.CurrentVideoReader = VideoReader(char(app.CurrentFile));
                app.CurrentVideoFile = app.CurrentFile;
                ok = true;
            catch exc
                uialert(app.UIFigure, exc.message, 'Video Read Error');
            end
        end

        function framePool = getFramePool(~, videoReader, n)
            totalFrames = round(videoReader.Duration * videoReader.FrameRate);
            n = min(n, totalFrames);
            frameIdx = sort(randperm(totalFrames, n));
            framePool = zeros(videoReader.Height, videoReader.Width, n, sprintf('uint%d', videoReader.BitsPerPixel));
            for idx = 1:n
                try
                    framePool(:, :, idx) = read(videoReader, frameIdx(idx));
                catch
                end
            end
        end

        function idx = bestRegion(~, stats, refX, refY, targetArea)
            n = numel(stats);
            if n == 1
                idx = 1;
                return;
            end
            centroids = vertcat(stats.Centroid);
            dists = sqrt((centroids(:, 1) - refX).^2 + (centroids(:, 2) - refY).^2);
            areas = [stats.Area]';
            intensities = [stats.MeanIntensity]';
            [~, order] = sort(dists, 'ascend');
            rankPos = zeros(n, 1);
            rankPos(order) = 1:n;
            [~, order] = sort(abs(areas - targetArea), 'ascend');
            rankArea = zeros(n, 1);
            rankArea(order) = 1:n;
            [~, order] = sort(intensities, 'descend');
            rankIntensity = zeros(n, 1);
            rankIntensity(order) = 1:n;
            [~, idx] = min(rankPos + rankArea + rankIntensity);
        end
    end
end
%}