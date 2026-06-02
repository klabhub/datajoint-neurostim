classdef PupilTracker < handle
    % sbx.PupilTracker  Interactive pupil tracker for Scanbox eye movies.
    %
    % QUICK START
    %   pt = sbx.PupilTracker(experiment)   % create object and queue files
    %   pt.initialize()                     % interactive parameter setup
    %   pt.track()                          % run automated tracking
    %
    % CONSTRUCTION
    %   pt = sbx.PupilTracker(experiment)
    %   pt = sbx.PupilTracker(experiment, overwrite=true)
    %   pt = sbx.PupilTracker(experiment, initialize=true)
    %
    %   'experiment' can be:
    %     - an ns.Experiment table row
    %     - a folder/subject path (e.g. "data\2024\02\23\252")
    %     - the full path to an existing _pupil.json parameter file
    %   The constructor locates all associated _eye.mj2 files. If a
    %   _pupil.json already exists its parameters are loaded automatically.
    %   Pass overwrite=true to force re-initialization.
    %
    %  When scanning folders for files, the pupiltracker will ask whether
    %  existing configuration json files (the output of initialize) should
    %  be overwritten, skipped, or read. The same options are given for
    % existing tsv files (the output of track).
    %
    % WORKFLOW
    %   Step 1 — initialize()
    %     Samples MaxSampleFrames random frames, shows their pixel-average,
    %     and asks you to draw a blue ellipse around the eye region.
    %
    %     GOAL: draw an ellipse where the pupil might appear. Try to mark
    %     the eye ball, and avoid the eye lid because that is often quite
    %     bright  and could be mistaken for the pupil. If the image is
    %     unclear, press the button to read a new random set of 100 frames.
    %     Sometimes pressing the "new average" button a few times is needed
    %     to get a good idea of the outline of the eye ball. Double click
    %     or press enter to continue.
    %
    %     To skip a file, close the figure (no json file with
    %     initialization parameters will be generated for that movie).
    %
    %     The figure now shows 6 (NrPupilImages) individual frames (brightness-ranked,
    %     blinks discarded). If the blue eye outline looks wrong (e.g., it
    %    overlaps with the eye lid), press restart to redo this file.
    %
    %     GOAL: draw a red ellipse to capture the (bright) pupil as
    %     closely as you can. Confirm each one with a double click or by pressing
    %     enter. If no pupil is visible, then draw the red "pupil" ellipse
    %     outside the blue eye ellipse; the subsequent code will ignore that image.
    %
    %     Once all have been confirmed, parameters are written to a JSON file next to
    %     the .mj2 movie so they can be reloaded without repeating this
    %     initialize() step.  Once this step is complete (for all experiments),
    %     the figure will show a three button interface to start the
    %     (automated) tracking.
    %
    %
    %   Step 2 — track()
    %     Reads every frame (or every frameStep-th frame), thresholds,
    %     morphologically cleans the binary image, and fits both a region-
    %     props ellipse and a bounding-box blob. Results are written to a
    %     tab-separated .tsv file next to the .mj2 movie.
    %
    %    To get an idea of how well the tracking works, click on Track
    %     with Preview. You can close the preview figure at any time;
    %     tracking will continue (and be much faster).
    %
    %
    % RESULTS TABLE COLUMNS
    %   Frame        – original video frame index
    %   X, Y         – pupil centroid (pixels)
    %   Area         – pupil area (pixels)
    %   MajorAxis    – ellipse major axis length (pixels)
    %   MinorAxis    – ellipse minor axis length (pixels)
    %   Eccentricity – ellipse eccentricity (0=circle, 1=line)
    %   Orientation  – ellipse orientation (degrees)
    %   BBox         – bounding box [x y w h] of the blob
    %   Threshold    – threshold used for this frame (reserved, currently NaN)
    %   FitQuality   – ratio of blob area to fitted ellipse area (1=perfect)
    %   EllipseIR    – ellipse mean intensity / mean eye intensity
    %   BlobIR       – blob mean intensity / mean eye intensity
    %
    % nsMeta INTEGRATION
    %  In nsMeta, pressing 'P' with a session selected will construct a PupilTracker object
    %  for that session's eye movies and launch the interactive initialization routine.
    %
    % DATAJOINT INTEGRATION
    %   populate(ns.C,'ctag="pupil"') will:
    %     * raise an error (with link) if no JSON exists yet
    %     * run track() if a JSON exists but no TSV
    %     * load the TSV into the CChannel table if both files exist
    %
    %   Example CParm setup:
    %     ptParms = struct('method','PupilTracker', ...
    %                      'variables',{{"X","Y","Area","FitQuality", ...
    %                                    "Threshold","IntensityRatio"}});
    %     eyePrep = struct('ctag','pupil','extension','.mj2', ...
    %                      'include','%_eye.mj2','fun','sbx.readMovie', ...
    %                      'description','Use sbx.PupilTracker', ...
    %                      'parms',ptParms);
    %     insertIfNew(ns.CParm, eyePrep);
    %
    % See also  sbx.PupilTracker.initialize, sbx.PupilTracker.track,
    %           sbx.PupilTracker.read, sbx.PupilTracker.plot
    %
    % BK - Mar 2026
    properties (Constant)
        FigureTag = 'sbx.PupilTracker.Figure'  % Unique name of the reused figure
    end

    properties
        Parameters  % Map (filename -> parameter row table) for tracking, determined with the interactive initialization routine.
        Results     % Map (filename -> results table) with tracking results for each video
        Figure      % Handle to the GUI figure

        MaxSampleFrames = 100  % Frames used in the interactive initialization routine
        MinRadius       = 2   % Pixels used with imopen to clean small elements.
        MinAreaFrac     = 0.05  % Fraction of the median pupil area drawn by the user.
        NrPupilImages   = 6;   % Number of frames shown in the interactive initialization routine for pupil selection (ranked by brightness within the eye)
    end

    methods
        function obj = PupilTracker(experiment, pv)
            % PUPIL Construct an instance of the pupil tracking class
            %   obj = sbx.PupilTracker(experiment)
            arguments
                experiment (1,1)
                pv.endsWith (1,1) string = "_eye.mj2"
                pv.overwrite (1,1) logical = false
                pv.autoskip (1,1) logical  =false
                pv.initialize (1,1) logical = false  % Open figure to initialize
                pv.modal (1,1) logical = false % Block execution until initialization is complete
            end
            if isa(experiment,"ns.Experiment")
                % Extract files from datajoint experiment table
                files = experiment * (ns.File & sprintf('filename LIKE "%%%s"',pv.endsWith));
                movieFiles = fullfile(folder(experiment & proj(files)), fetchn(files, 'filename'));
            elseif endsWith(experiment,"_pupil.json")
                % experiment is a json parameter file previously produced by
                % PupilTracker. Read params and run track non-interactively.
                jsonPath = char(experiment);
                paramStruct = jsondecode(fileread(jsonPath));
                % Derive the .mj2 path from the JSON filename and the stored FileName field
                [jsonDir, ~, ~] = fileparts(jsonPath);
                mj2File = fullfile(jsonDir, strrep(paramStruct.FileName, '_pupil.json', '.mj2'));
                movieFiles = string(mj2File);

                obj.Parameters = containers.Map('KeyType', 'char', 'ValueType', 'any');
                obj.Results    = containers.Map('KeyType', 'char', 'ValueType', 'any');
                paramStruct = rmfield(paramStruct, 'FileName');
                obj.Parameters(char(movieFiles)) = struct2table(paramStruct);
                fprintf('Loaded parameters from: %s\n', jsonPath);
                track(obj,visualize=false);
                return;
            elseif endsWith(experiment,"_eye.mj2")
                movieFiles = experiment;
            else
                % experiment is a string with a session name (from nsMeta)
                % e,g, \root\2024\02\23\252
                files = dir(experiment + "*\**\*" + pv.endsWith);
                movieFiles = string(fullfile({files.folder},{files.name}))';
            end

            % Pre-populate parameters and results maps (keyed by full file path)

            numFiles = numel(movieFiles);
            jsonFiles = strrep(movieFiles, '.mj2', '_pupil.json'); % Files with initialization parameters

            %% Check for existing parameter JSON files
            skip = false(numFiles,1);
            readJson = false(numFiles,1);
            if pv.overwrite
                overwrite = true(numFiles,1);
            else
                overwrite = false(numFiles,1);
                for i = 1:numFiles
                    if exist(char(jsonFiles(i)), 'file')
                        if pv.autoskip
                            skip(i) =true;
                        else
                            response = questdlg( ...
                                sprintf('Parameter file already exists:\n%s\n\nWhat would you like to do?', jsonFiles(i)), ...
                                'Existing Parameters', ...
                                'Redo', 'Skip', 'Read', 'Read');
                            if isempty(response), response = 'Skip'; end
                            switch response
                                case 'Overwrite'
                                    overwrite(i) = true;
                                case 'Skip'
                                    skip(i) = true;
                                case 'Read'
                                    readJson(i) = true;
                            end
                        end
                    end
                end
            end

            jsonFiles(skip) = [];
            movieFiles(skip) = [];
            overwrite(skip) = [];
            readJson(skip) = [];
            numFiles = numel(movieFiles);
            readTsv = false(numFiles,1);

            %% Check for existing output (TSV) files
            outputFiles = strrep(movieFiles, '.mj2', '_pupil.tsv'); % Files with tracking results
            if ~pv.overwrite
                skip  = false(numFiles,1);

                for i = 1:numFiles
                    if exist(char(outputFiles(i)), 'file')
                        if pv.autoskip
                            skip(i) = true;
                        else
                            response = questdlg( ...
                                sprintf('Output file already exists:\n%s\n\nWhat would you like to do?', outputFiles(i)), ...
                                'Existing Output', ...
                                'Overwrite', 'Skip', 'Read', 'Skip');
                            if isempty(response)
                                % Dialog closed — treat as cancel
                                error('Execution cancelled by user.');
                            end
                            switch response
                                case 'Overwrite'
                                    overwrite(i) = true;
                                case 'Skip'
                                    skip(i) = true;
                                case 'Read'
                                    readTsv(i) = true;
                            end
                        end
                    end
                end
            end

            jsonFiles(skip) = [];
            movieFiles(skip) = [];
            outputFiles(skip) = [];
            overwrite(skip) = [];
            readJson(skip) = [];
            readTsv(skip) = [];


            emptyRow = table(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
                'VariableNames', {'Threshold', 'Floor', 'StartX', 'StartY', 'Area', ...
                'EyeX', 'EyeY', 'EyeMajor', 'EyeMinor', 'EyeOrientation'});
            obj.Parameters = containers.Map('KeyType', 'char', 'ValueType', 'any');
            obj.Results    = containers.Map('KeyType', 'char', 'ValueType', 'any');
            numFiles = numel(movieFiles);
            for i = 1:numFiles
                obj.Parameters(char(movieFiles(i))) = emptyRow;
            end
            % Read json
            for i= 1:numFiles
                if readJson(i)
                    try
                        paramStruct = jsondecode(fileread(char(jsonFiles(i))));
                        paramStruct = rmfield(paramStruct, 'FileName');
                        obj.Parameters(char(movieFiles(i))) = struct2table(paramStruct);
                        fprintf('  -> Loaded: %s\n', jsonFiles(i));
                    catch
                        warning('Could not read %s. Parameters left empty.', jsonFiles(i));
                    end
                end
                if readTsv(i)
                    obj.Results(char(movieFiles(i))) = readtable(char(outputFiles(i)), 'FileType', 'text', 'Delimiter', '\t');
                    fprintf('  -> Loaded: %s\n', outputFiles(i));
                end
            end

            if pv.initialize
                initialize(obj,modal=pv.modal);
            end
        end
        function disp(obj)
            % DISP - Display a summary of the object.
            paramFiles = keys(obj.Parameters);
            fprintf('sbx.PupilTracker with %d file(s) queued for processing:\n', numel(paramFiles));
            for i = 1:numel(paramFiles)
                fprintf('  [%d] %s\n', i, paramFiles{i});
            end
        end
        function delete(obj)
            % DELETE - Close the figure when the object is deleted.
            if ~isempty(obj.Figure)  && isgraphics(obj.Figure)
                close(obj.Figure);
            end
        end
        function initialize(obj,pv)
            % INITIALIZE  Interactively set tracking parameters for all queued videos.
            %
            %   obj.initialize()
            %   obj.initialize(overwrite=true)
            %
            % For each video the routine proceeds in two phases:
            %
            %   Phase 1 — Eye boundary
            %     Displays the pixel-average of MaxSampleFrames random frames.
            %     Draw a blue ellipse around the eye and double-click to confirm.
            %     "New Average" resamples and recomputes the average in-place.
            %     Close the window to skip that file.
            %
            %   Phase 2 — Pupil selection
            %     The same frame pool is blink-filtered (frames whose mean
            %     eye-region brightness falls below the 5th percentile of the
            %     average frame are discarded) and then NrPupilImages frames
            %     are chosen at equal intervals across the remaining brightness
            %     rank. Draw a red ellipse over the pupil in each tile.
            %     Drawing outside the blue eye ellipse silently skips that tile.
            %
            % Parameters saved to JSON (reloaded automatically on next construction):
            %   Threshold      – detection threshold (midpoint of pupil / eye medians)
            %   Floor          – median eye-region intensity
            %   StartX, StartY – median pupil centroid across drawn ellipses
            %   Area           – median geometric area of drawn ellipses
            %   EyeX, EyeY     – eye ellipse centre
            %   EyeMajor, EyeMinor, EyeOrientation – eye ellipse dimensions
            %
            % Options:
            %   overwrite (logical, default false) — re-initialize files that
            %     already have parameters.
            %
            % After all files are processed an optional confirmation dialog
            % offers to run track() immediately (with or without preview).
            %
            % See also  sbx.PupilTracker, sbx.PupilTracker.track
            arguments
                obj
                pv.overwrite (1,1) logical = false
                pv.modal (1,1) logical = false
            end

            disp('Starting batch parameter initialization...');
            paramFiles = keys(obj.Parameters);
            for i = 1:numel(paramFiles)
                currentFile = paramFiles{i};
                jsonFileName = strrep(currentFile, '.mj2', '_pupil.json');

                % Skip files whose parameters were already loaded (e.g. from constructor)
                p = obj.Parameters(currentFile);
                if p.Threshold ~= 0  && ~ pv.overwrite
                    fprintf('  [%d/%d] Parameters already set, skipping: %s\n', i, numel(paramFiles), currentFile);
                    continue;
                end

                fprintf('Processing %d of %d: %s\n', i, numel(paramFiles), currentFile);
                try
                    v = VideoReader(char(currentFile));
                catch
                    warning('Could not read file. Skipping...');
                    continue;
                end

                [~, fname, ext] = fileparts(currentFile);
                obj.Figure = findobj('Type', 'figure', 'Tag', obj.FigureTag);
              
                if isempty(obj.Figure)
                    obj.Figure = figure('Name', [fname ext], 'Tag', obj.FigureTag, 'units', 'normalized', 'position', [0.05 0.05 0.9 0.9]);
                else
                    figure(obj.Figure); % Bring to front
                    clf(obj.Figure);
                    obj.Figure.Name = [fname ext];
                end
                if pv.modal
                    set(obj.Figure, 'WindowStyle', 'modal');
                end

                %% ---- Sample frames once (reused for both eye and pupil selection) ----
                skippedFile = false;
                restartFile = true;
                while restartFile
                    restartFile = false;
                    setappdata(obj.Figure, 'restartFile', false);

                    framePool = sbx.PupilTracker.getFramePool(v, obj.MaxSampleFrames);
                    brightnessPerFrame = squeeze(mean(framePool,[1 2]));
                    keep = brightnessPerFrame > prctile(brightnessPerFrame,50);
                    eyeSelFrame = uint8(mean(double(framePool(:,:,keep)), 3));

                    %% ---- Phase 1: Eye selection on averaged frame -------------------------
                    clf(obj.Figure);
                    ax_eye = axes(obj.Figure, 'Position', [0 0.07 1 0.93]);
                    hImg = imshow(eyeSelFrame, 'Parent', ax_eye, 'InitialMagnification', 'fit');
                    title(ax_eye, 'Draw eye boundary (blue). Double-click to confirm. Close window to skip.', ...
                        'Color', 'black', 'FontSize', 10);
                    uicontrol(obj.Figure, 'Style', 'pushbutton', 'String', 'New Average', ...
                        'Units', 'normalized', 'Position', [0.01 0.01 0.15 0.05], ...
                        'Callback', @(~,~) sbx.PupilTracker.resampleAverageFrame(v, ax_eye, obj.MaxSampleFrames));
                    setappdata(obj.Figure, 'skipFile', false);
                    uicontrol(obj.Figure, 'Style', 'pushbutton', 'String', 'Skip File', ...
                        'Units', 'normalized', 'Position', [0.17 0.01 0.15 0.05], ...
                        'BackgroundColor', [0.6 0.6 0.1], 'ForegroundColor', 'white', ...
                        'FontSize', 10, 'FontWeight', 'bold', ...
                        'Callback', @(~,~) sbx.PupilTracker.resumeWithFlag(obj.Figure, 'skipFile', true));

                    eyeRoi = drawellipse(ax_eye, 'Color', 'b', 'FaceAlpha', 0.15);
                    lEye = addlistener(eyeRoi, 'ROIClicked', @(~,evt) sbx.PupilTracker.roiConfirmIfDouble(obj.Figure, evt));
                    uiwait(obj.Figure);
                    delete(lEye);
                    if ~isempty(obj.Figure) && ~isgraphics(obj.Figure)
                        fprintf('  Skipped (figure closed): %s\n', currentFile);
                        remove(obj.Parameters, currentFile);
                        skippedFile = true; break;
                    end
                    if getappdata(obj.Figure, 'skipFile')
                        fprintf('  Skipped (button): %s\n', currentFile);
                        remove(obj.Parameters, currentFile);
                        skippedFile = true; break;
                    end
                    if getappdata(obj.Figure, 'restartFile'), restartFile = true; continue; end
                    eyeRoi.InteractionsAllowed = 'none';

                    % Capture eye parameters before clearing the figure
                    eyeCenter   = eyeRoi.Center;
                    eyeSemiAxes = eyeRoi.SemiAxes;
                    eyeAngle    = eyeRoi.RotationAngle;

                    % Build eye mask geometrically (matches track() sign convention)
                    [Xg, Yg]  = meshgrid(1:v.Width, 1:v.Height);
                    a_eye     = eyeSemiAxes(1);  b_eye = eyeSemiAxes(2);
                    theta_eye = -deg2rad(eyeAngle);
                    eyeMask   = (((Xg - eyeCenter(1)) * cos(theta_eye) + (Yg - eyeCenter(2)) * sin(theta_eye)).^2 / a_eye^2 + ...
                        ((Xg - eyeCenter(1)) * sin(theta_eye) - (Yg - eyeCenter(2)) * cos(theta_eye)).^2 / b_eye^2) <= 1;



                    %% ---- Phase 2: rank frames by brightness inside the eye ---------------
                    fprintf('  Ranking %d frames by eye-region brightness...\n', size(framePool, 3));
                    nPool     = size(framePool, 3);

                    % Discard frames whose mean eye-region intensity is below the 5th percentile of the
                    % averaged eye frame.  These are probably blinks
                    eyeSelFrame = hImg(1).CData;
                    eyePixelsRef  = double(eyeSelFrame(eyeMask));
                    eyeFloorRef   = prctile(eyePixelsRef,5);
                    eyeBrightness = zeros(nPool, 1);
                    keepFrame     = true(nPool, 1);
                    for f = 1:nPool
                        frameSlice       = double(framePool(:,:,f));
                        eyeBrightness(f) = mean(frameSlice(eyeMask));
                        if eyeBrightness(f) < eyeFloorRef
                            keepFrame(f) = false;
                        end
                    end
                    nDiscarded = sum(~keepFrame);
                    if nDiscarded > 0
                        fprintf('  Discarded %d frame(s) as likely blinks (eye mean < 5th percentile).\n', nDiscarded);
                    end
                    % Rerank only the kept frames; fall back to all frames if too few survive
                    keptIdx = find(keepFrame);
                    if numel(keptIdx) < obj.NrPupilImages
                        warning('Too few frames passed the blink filter (%d); using all frames.', numel(keptIdx));
                        keptIdx = 1:nPool;
                    end
                    [~, rankAmongKept] = sort(eyeBrightness(keptIdx));   % ascending within kept set
                    rankOrder          = keptIdx(rankAmongKept);

                    % Pick equally distributed over the rank range

                    nKept    = numel(rankOrder);
                    quantPos = round(linspace(1, nKept, obj.NrPupilImages));
                    tileIdx  = rankOrder(quantPos);
                    framePool = framePool(:,:,tileIdx);
                    % Eye outline coordinates for overlay (same rotation as track())
                    phi_eye = linspace(0, 2*pi, 100);
                    alpha   = deg2rad(eyeAngle);
                    R_eye   = [cos(alpha) sin(alpha); -sin(alpha) cos(alpha)];
                    xy_eye  = R_eye * [a_eye * cos(phi_eye); b_eye * sin(phi_eye)];

                    %% Show tiled frames with eye outline
                    clf(obj.Figure);
                    tLayout = tiledlayout(obj.Figure, "Flow", 'TileSpacing', 'compact', 'Padding', 'compact');
                    sampleAxes   = cell(obj.NrPupilImages, 1);
                    for t = 1:obj.NrPupilImages
                        sampleAxes{t}   = nexttile(tLayout);
                        imshow(framePool(:,:,t), 'Parent', sampleAxes{t}, 'InitialMagnification', 'fit');
                        hold(sampleAxes{t}, 'on');
                        plot(sampleAxes{t}, xy_eye(1,:) + eyeCenter(1), xy_eye(2,:) + eyeCenter(2), ...
                            'b-', 'LineWidth', 1.5);
                        hold(sampleAxes{t}, 'off');
                        title(sampleAxes{t}, ...
                            sprintf('Rank %d/%d — draw pupil (red). Double-click to confirm.', quantPos(t), nKept), ...
                            'Color', 'black', 'FontSize', 7);
                    end
                    drawnow;

                    % Draw pupil ellipses one by one (Restart/Skip buttons above tile grid)
                    uicontrol(obj.Figure, 'Style', 'pushbutton', 'String', 'Restart', ...
                        'Units', 'normalized', 'Position', [0.84 0.01 0.15 0.05], ...
                        'BackgroundColor', [0.75 0.25 0.1], 'ForegroundColor', 'white', ...
                        'FontSize', 10, 'FontWeight', 'bold', ...
                        'Callback', @(~,~) sbx.PupilTracker.restartInitFile(obj.Figure));
                    setappdata(obj.Figure, 'skipFile', false);
                    uicontrol(obj.Figure, 'Style', 'pushbutton', 'String', 'Skip File', ...
                        'Units', 'normalized', 'Position', [0.68 0.01 0.15 0.05], ...
                        'BackgroundColor', [0.6 0.6 0.1], 'ForegroundColor', 'white', ...
                        'FontSize', 10, 'FontWeight', 'bold', ...
                        'Callback', @(~,~) sbx.PupilTracker.resumeWithFlag(obj.Figure, 'skipFile', true));
                    skipFile  = false;
                    pupilRois = cell(obj.NrPupilImages, 1);
                    for t = 1:obj.NrPupilImages
                        pupilRois{t} = drawellipse(sampleAxes{t}, 'Color', 'r', 'FaceAlpha', 0.15);
                        % Check flags immediately after drawellipse's blocking draw phase.
                        % If Restart/Skip was clicked while the user was drawing, uiresume could not
                        % interrupt drawellipse's internal event loop — but setappdata still ran.
                        % Catching the flag here avoids proceeding to uiwait unnecessarily.
                        if ~isgraphics(obj.Figure), skipFile = true; break; end
                        if getappdata(obj.Figure, 'restartFile'), restartFile = true; break; end
                        if getappdata(obj.Figure, 'skipFile'),    skipFile    = true; break; end
                        lPupil = addlistener(pupilRois{t}, 'ROIClicked', @(~,evt) sbx.PupilTracker.roiConfirmIfDouble(obj.Figure, evt));
                        uiwait(obj.Figure);
                        delete(lPupil);
                        if ~isgraphics(obj.Figure), skipFile = true; break; end
                        if getappdata(obj.Figure, 'restartFile'), restartFile = true; break; end
                        if getappdata(obj.Figure, 'skipFile'),    skipFile    = true; break; end
                        pupilRois{t}.InteractionsAllowed = 'none';
                        title(sampleAxes{t}, sprintf('Rank %d/%d', quantPos(t), nPool), 'Color', 'w', 'FontSize', 7);
                    end
                    if restartFile, continue; end  % while restartFile — restart from Phase 1
                    if skipFile || (~isempty(obj.Figure)  &&  ~isgraphics(obj.Figure))
                        fprintf('  Skipped (figure closed): %s\n', currentFile);
                        remove(obj.Parameters, currentFile);
                        skippedFile = true; break;
                    end

                    %% --- Extract and Calculate Pupil Parameters ---
                    % Collect pupil pixels only from ellipses whose centroid is inside the eye mask
                    % This allows users to skip one frame by drawing outside
                    % the blue ellipse
                    allPupilPixels = [];
                    allEyePixels = [];
                    validCentroids = zeros(0, 2);
                    validAreas     = zeros(0, 1);
                    for t = 1:obj.NrPupilImages
                        if  ~isgraphics(pupilRois{t}), continue; end
                        cx_p = pupilRois{t}.Center(1);  cy_p = pupilRois{t}.Center(2);
                        ci = round(cy_p);  cj = round(cx_p);
                        if ci < 1 || ci > v.Height || cj < 1 || cj > v.Width, continue; end
                        if ~eyeMask(ci, cj), continue; end   % centroid must lie inside the eye
                        maskT  = createMask(pupilRois{t});
                        thisFrame = framePool(:,:,t);
                        thisEyePixels = double(thisFrame(eyeMask & ~maskT));
                        thisPupilPixels = double(thisFrame(maskT));
                        allPupilPixels             = [allPupilPixels; thisPupilPixels];           %#ok<AGROW>
                        allEyePixels             = [allEyePixels             ; thisEyePixels];           %#ok<AGROW>
                        validCentroids(end+1, :)   = [cx_p, cy_p];                      %#ok<AGROW>
                        validAreas(end+1)           = pi * prod(pupilRois{t}.SemiAxes);  %#ok<AGROW>
                    end

                    if isempty(allPupilPixels)
                        warning('No valid pupil ellipses drawn inside the eye boundary for %s. Restarting...', currentFile);
                        restartFile = true; continue;  % while restartFile
                    end

                    % Store all parameters as a single row in the map
                    obj.Parameters(currentFile) = table( ...
                        (median(allPupilPixels)+median(allEyePixels))/2, ...            % Threshold: midpoint of pupil/eye medians
                        median(allEyePixels), ...                                      % EyeMedian: median eye-region intensity (threshold floor)
                        median(validCentroids(:, 1)), median(validCentroids(:, 2)), ...  % StartX, StartY: mean centroid of valid pupil ellipses
                        median(validAreas), ...                                      % TargetArea: mean geometric area of drawn pupil ellipses
                        eyeCenter(1), eyeCenter(2), ...                     % EyeCenterX, EyeCenterY
                        eyeSemiAxes(1) * 2, eyeSemiAxes(2) * 2, eyeAngle, ... % EyeMajorAxis, EyeMinorAxis, EyeOrientation
                        'VariableNames', {'Threshold', 'Floor','StartX', 'StartY', 'Area',...
                        'EyeX', 'EyeY', 'EyeMajor', 'EyeMinor', 'EyeOrientation'});

                    % Write parameters to a JSON file (store only filename, not full path)
                    paramStruct = table2struct(obj.Parameters(currentFile));
                    paramStruct.FileName = [fname ext];
                    jsonTxt = jsonencode(paramStruct, "PrettyPrint", true);
                    fid = fopen(jsonFileName, 'w');
                    if fid == -1
                        warning('Could not write parameters to %s', jsonFileName);
                    else
                        fwrite(fid, jsonTxt, 'char');
                        fclose(fid);
                        fprintf('  -> Saved parameters to: %s\n', jsonFileName);
                    end

                    %% Per-file track dialog
                    clf(obj.Figure);
                    setappdata(obj.Figure, 'runTrack', false);
                    setappdata(obj.Figure, 'runTrackPreview', false);
                    setappdata(obj.Figure, 'nextFile', false);
                    obj.Figure.Name = sprintf('[%d/%d] %s initialized', i, numel(paramFiles), [fname ext]);
                    uicontrol(obj.Figure, 'Style', 'text', ...
                        'String', sprintf('[%d/%d]  %s  — saved.', i, numel(paramFiles), [fname ext]), ...
                        'Units', 'normalized', 'Position', [0.1 0.64 0.8 0.08], ...
                        'FontSize', 11, 'HorizontalAlignment', 'center');
                    uicontrol(obj.Figure, 'Style', 'pushbutton', 'String', 'Track with Preview', ...
                        'Units', 'normalized', 'Position', [0.35 0.53 0.30 0.09], ...
                        'BackgroundColor', [0.2 0.5 0.8], 'ForegroundColor', 'white', ...
                        'FontSize', 12, 'FontWeight', 'bold', ...
                        'Callback', @(~,~) sbx.PupilTracker.resumeWithFlag(obj.Figure, 'runTrackPreview', true));
                    uicontrol(obj.Figure, 'Style', 'pushbutton', 'String', 'Track', ...
                        'Units', 'normalized', 'Position', [0.35 0.42 0.30 0.09], ...
                        'BackgroundColor', [0.2 0.6 0.2], 'ForegroundColor', 'white', ...
                        'FontSize', 12, 'FontWeight', 'bold', ...
                        'Callback', @(~,~) sbx.PupilTracker.resumeWithFlag(obj.Figure, 'runTrack', true));
                    if i < numel(paramFiles)
                        nextLabel = sprintf('Next File  (%d / %d)', i+1, numel(paramFiles));
                    else
                        nextLabel = 'Done';
                    end
                    uicontrol(obj.Figure, 'Style', 'pushbutton', 'String', nextLabel, ...
                        'Units', 'normalized', 'Position', [0.35 0.31 0.30 0.08], ...
                        'FontSize', 11, ...
                        'Callback', @(~,~) sbx.PupilTracker.resumeWithFlag(obj.Figure, 'nextFile', true));
                    uiwait(obj.Figure);

                    % Track this file only if requested
                    doPreview = isgraphics(obj.Figure) && isappdata(obj.Figure, 'runTrackPreview') && getappdata(obj.Figure, 'runTrackPreview');
                    doTrack   = isgraphics(obj.Figure) && isappdata(obj.Figure, 'runTrack')        && getappdata(obj.Figure, 'runTrack');
                    if doPreview || doTrack
                        savedParams    = obj.Parameters;
                        obj.Parameters = containers.Map({currentFile}, {savedParams(currentFile)});
                        if doPreview
                            obj.track(visualize=true, frameStep=10, save=false);
                            obj.Parameters = savedParams;
                            % Preview window was closed by the user — they are done; skip
                            % remaining files and the completion menu.
                            if isempty(obj.Figure) || ~isgraphics(obj.Figure)
                                return;
                            end
                        else
                            obj.track(visualize=false);
                            obj.Parameters = savedParams;
                        end
                    end
                end  % while restartFile
                if skippedFile, continue; end  % for i — skip to next file
            end
            disp('Finished collecting parameters. Run track() to generate a pupil tracking tsv file.');
            % Only show the track menu if there are files with parameters ready to track.
            if obj.Parameters.Count == 0
                if ~isempty(obj.Figure) && isgraphics(obj.Figure), close(obj.Figure); end
                return;
            end
            % Ensure a figure exists for the completion dialog even when all files were
            % already initialized and skipped (obj.Figure was never created in the loop).
            if isempty(obj.Figure) || ~isgraphics(obj.Figure)
                obj.Figure = figure('Name', 'Initialization complete', 'Tag', obj.FigureTag, ...
                    'Units', 'normalized', 'Position', [0.05 0.05 0.9 0.9]);
                if pv.modal
                    set(obj.Figure, 'WindowStyle', 'modal');
                end
            end
            if isgraphics(obj.Figure)
                setappdata(obj.Figure, 'runTrack', false);
                setappdata(obj.Figure, 'runTrackPreview', false);
                obj.Figure.Name = 'Initialization complete';
                uicontrol(obj.Figure, 'Style', 'text', ...
                    'String', 'All files initialized.', ...
                    'Units', 'normalized', 'Position', [0.1 0.62 0.8 0.08], ...
                    'FontSize', 13, 'HorizontalAlignment', 'center');
                uicontrol(obj.Figure, 'Style', 'pushbutton', 'String', 'Track with Preview', ...
                    'Units', 'normalized', 'Position', [0.35 0.51 0.30 0.09], ...
                    'BackgroundColor', [0.2 0.5 0.8], 'ForegroundColor', 'white', ...
                    'FontSize', 12, 'FontWeight', 'bold', ...
                    'Callback', @(~,~) sbx.PupilTracker.resumeWithFlag(obj.Figure, 'runTrackPreview', true));
                uicontrol(obj.Figure, 'Style', 'pushbutton', 'String', 'Track', ...
                    'Units', 'normalized', 'Position', [0.35 0.40 0.30 0.09], ...
                    'BackgroundColor', [0.2 0.6 0.2], 'ForegroundColor', 'white', ...
                    'FontSize', 12, 'FontWeight', 'bold', ...
                    'Callback', @(~,~) sbx.PupilTracker.resumeWithFlag(obj.Figure, 'runTrack', true));
                uicontrol(obj.Figure, 'Style', 'pushbutton', 'String', 'Close', ...
                    'Units', 'normalized', 'Position', [0.35 0.29 0.30 0.08], ...
                    'FontSize', 11, ...
                    'Callback', @(~,~) close(obj.Figure));
                uiwait(obj.Figure);
            end
            if ~isempty(obj.Figure) && isgraphics(obj.Figure) && ...
                    isappdata(obj.Figure, 'runTrackPreview') && getappdata(obj.Figure, 'runTrackPreview')
                obj.track(visualize=true, frameStep=10, save=false);
            elseif ~isempty(obj.Figure)  && isgraphics(obj.Figure) && ...
                    isappdata(obj.Figure, 'runTrack') && getappdata(obj.Figure, 'runTrack')
                close(obj.Figure);
                obj.Figure = [];
                obj.track();
            end
        end

        function track(obj, pv)
            % TRACK  Run automated pupil tracking on all queued videos.
            %
            %   obj.track()
            %   obj.track(visualize=true)
            %   obj.track(maxFrames=500)
            %   obj.track(frameStep=3)
            %
            % For each video the method:
            %   1. Reads frames indexed by frameIndices = 1:frameStep:NumFrames
            %   2. Applies the threshold from initialize() within the eye mask
            %   3. Cleans the binary image with imopen + imfill + bwareaopen
            %   4. Calls regionprops to fit an ellipse and a bounding-box blob
            %   5. Selects the best region via composite rank (centroid distance,
            %      area deviation, mean intensity) using bestRegion()
            %   6. Stores all results as a tab-separated .tsv file next to the
            %      .mj2 movie and in obj.Results
            %
            % Options:
            %   visualize (logical, default false) — show each processed frame
            %     with the fitted ellipse and bounding box overlaid.
            %   maxFrames (double, default inf) — stop after this many processed
            %     frames per video (useful for quick sanity checks).
            %   frameStep (double, default 1) — process every n-th frame;
            %     e.g. frameStep=5 processes frames 1, 6, 11, … and runs ~5x
            %     faster at the cost of temporal resolution.
            %
            % Output columns — see help(sbx.PupilTracker) for the full list.
            %
            % See also  sbx.PupilTracker, sbx.PupilTracker.initialize,
            %           sbx.PupilTracker.read, sbx.PupilTracker.plot
            arguments
                obj
                pv.visualize (1,1) logical = false  % Set to true to show the frames and the detected regions
                pv.maxFrames (1,1) double = inf  % For quick debugging; process only up to maxFrames
                pv.frameStep (1,1) double = 1    % Step size between processed frames (>1 fast-forwards)
                pv.save (1,1) logical = true     % Set to false to skip writing results to disk
            end

            if obj.Parameters.Count == 0
                disp('All files were skipped. Nothing to process.');
                return;
            end

            % Check if parameters have been initialized (Threshold is 0 for all entries
            % when initialize() has not yet been run).
            paramVals = values(obj.Parameters);
            if all(cellfun(@(p) p.Threshold == 0, paramVals))
                fprintf('No parameters found. Running initialize() first...\n');
                obj.initialize();
            end

            videoFiles = keys(obj.Parameters);
            numVideos = numel(videoFiles);
            fprintf('Starting ellipse and blob tracking for %d videos...\n', numVideos);

            if pv.visualize
                obj.Figure = findobj('Type', 'figure', 'Tag', obj.FigureTag);
                if isempty(obj.Figure)
                    obj.Figure = figure('Name', 'Pupil Tracker Viewer', 'Tag', obj.FigureTag, 'Units', 'Normalized', 'Position', [0.25 0.25 0.5 0.5]);
                else
                    figure(obj.Figure); % Bring to front
                    clf(obj.Figure);
                    obj.Figure.Name = 'Pupil Tracker Viewer';
                end
                % Create preview axes leaving space at the bottom for the scrubber and top for the title
                axPreview = axes(obj.Figure, 'Position', [0 0.11 1 0.84], 'Tag', 'previewAxes');
                setappdata(obj.Figure, 'previewAxesHandle', axPreview);
                setappdata(obj.Figure, 'previewMode', 'tracking');
                setappdata(obj.Figure, 'previewLastRenderedFrame', NaN);
                setappdata(obj.Figure, 'previewIgnoreSeekCallback', false);
                % Position scrubber — use a Java JScrollBar so the thumb width is
                % proportional to ~5% of the total range (easier to grab)
                setappdata(obj.Figure, 'previewSeekFrame', 0);
                jSB = javax.swing.JScrollBar(javax.swing.JScrollBar.HORIZONTAL);
                jSB.setValues(1, 1, 1, 2);  % (value, extent, min, max) — updated per video
                ws = warning('off', 'MATLAB:javacomponent:FunctionToBeRemoved');
                [~, hScrubberCont] = javacomponent(jSB, [], obj.Figure);
                warning(ws);
                hScrubberCont.Units = 'normalized';
                hScrubberCont.Position = [0.01 0.01 0.88 0.09];
                hScrubberCont.Tag = 'previewScrubber';
                setappdata(obj.Figure, 'previewJScrollBar', jSB);
                set(handle(jSB, 'CallbackProperties'), 'AdjustmentValueChangedCallback', ...
                    @(~,~) sbx.PupilTracker.seekToFrame(obj.Figure, jSB));
                uicontrol(obj.Figure, 'Style', 'text', 'String', '0 / 0', ...
                    'Units', 'normalized', 'Position', [0.90 0.01 0.09 0.10], ...
                    'Tag', 'previewFrameLabel', 'HorizontalAlignment', 'center', 'FontSize', 9);
            end

            trackStart = tic; %#ok<NASGU>
            for vIdx = 1:numVideos
                videoFile = videoFiles{vIdx};
                outputFile = strrep(videoFile, '.mj2', '_pupil.tsv');
                thisParameters = obj.Parameters(videoFile);

                se = strel('disk', obj.MinRadius);

                fprintf('\n[%d/%d] Processing: %s\n \t\t -> %s\n', vIdx, numVideos, videoFile, outputFile);

                try
                    v = VideoReader(char(videoFile));
                catch
                    warning('Could not read %s. Skipping...', videoFile);
                    continue;
                end

                lastFrame = min(v.NumFrames, floor(pv.maxFrames));
                videoStart = tic;

                % Update scrubber range for this video (thumb = ~5% of range)
                if pv.visualize && ~isempty(obj.Figure) && isgraphics(obj.Figure)
                    jSB = getappdata(obj.Figure, 'previewJScrollBar');
                    if ~isempty(jSB)
                        ext = max(1, round(0.05 * max(lastFrame - 1, 1)));
                        setappdata(obj.Figure, 'previewIgnoreSeekCallback', true);
                        jSB.setValues(1, ext, 1, lastFrame - 1 + ext);
                        setappdata(obj.Figure, 'previewIgnoreSeekCallback', false);
                    end
                    lbl = findobj(obj.Figure, 'Tag', 'previewFrameLabel');
                    if ~isempty(lbl), lbl(1).String = sprintf('1 / %d', lastFrame); end
                    setappdata(obj.Figure, 'previewSeekFrame', 0);
                    setappdata(obj.Figure, 'previewLastRenderedFrame', NaN);
                end

                % Core Arrays (pre-allocated for worst case; trimmed after loop)
                Frames = NaN(lastFrame, 1);
                X = NaN(lastFrame, 1);
                Y = NaN(lastFrame, 1);
                PupilArea = NaN(lastFrame, 1);
                ThresholdUsed = NaN(lastFrame, 1);

                % Ellipse-specific Arrays
                MajorAxis = NaN(lastFrame, 1);
                MinorAxis = NaN(lastFrame, 1);
                Eccentricity = NaN(lastFrame, 1);
                Orientation = NaN(lastFrame, 1);
                FitQuality = NaN(lastFrame, 1);
                EllipseIR = NaN(lastFrame, 1);
                BlobIR    = NaN(lastFrame, 1);

                % Blob-specific Arrays
                BoundingBox = NaN(lastFrame, 4);

                [Xgrid, Ygrid] = meshgrid(1:v.Width, 1:v.Height);

                % Create a static elliptical mask for the eye
                a = thisParameters.EyeMajor / 2;
                b = thisParameters.EyeMinor / 2;
                theta = -deg2rad(thisParameters.EyeOrientation); % Convert to radians and invert for calculation

                % Rotated ellipse equation
                eyeMask = (((Xgrid - thisParameters.EyeX) * cos(theta) + (Ygrid - thisParameters.EyeY) * sin(theta)).^2 / a^2 + ...
                    ((Xgrid - thisParameters.EyeX) * sin(theta) - (Ygrid - thisParameters.EyeY) * cos(theta)).^2 / b^2) <= 1;

                % Precompute eye outline for visualisation
                if pv.visualize
                    phi_eye_v = linspace(0, 2*pi, 50);
                    R_eye_v   = [cos(-theta) sin(-theta); -sin(-theta) cos(-theta)];
                    xy_eye_v  = R_eye_v * [a*cos(phi_eye_v); b*sin(phi_eye_v)];
                end

                progressStr = '';
                figClosed = false;
                iFrame = 0;
                frameNum = 1;

                while frameNum <= lastFrame
                    iFrame = min(iFrame + 1, lastFrame);  % cap to avoid overflow on repeated seeks
                    Frames(iFrame) = frameNum;
                    % Check whether the user dragged the scrubber to seek
                    if pv.visualize && ~isempty(obj.Figure) && isgraphics(obj.Figure)
                        seekTo = getappdata(obj.Figure, 'previewSeekFrame');
                        if seekTo > 0
                            frameNum = seekTo;
                            Frames(iFrame) = frameNum;
                            setappdata(obj.Figure, 'previewSeekFrame', 0);
                        end
                    end
                    %% Abort if figure was closed during preview
                    if pv.visualize && (isempty(obj.Figure) || ~isgraphics(obj.Figure))
                        fprintf('\n  Preview figure closed — stopping tracking for this file.\n');
                        figClosed = true;
                        break;
                    end
                    %% Progress indicator
                    if mod(iFrame, 25) == 0 || iFrame == 1
                        pctFrames  = 100 * frameNum / lastFrame;
                        pctOverall = 100 * (vIdx - 1 + frameNum / lastFrame) / numVideos;
                        elapsedVideo = toc(videoStart);
                        fpsEst = iFrame / max(elapsedVideo, 0.001);
                        etaSec = (lastFrame - frameNum) / max(pv.frameStep, 1) / fpsEst + ...
                            (numVideos - vIdx) * (lastFrame / fpsEst);
                        newMsg = sprintf('  [File %d/%d]  Frame: %d/%d (%.0f%%)  |  Overall: %.0f%%  |  ETA: %ds', ...
                            vIdx, numVideos, frameNum, lastFrame, pctFrames, pctOverall, round(etaSec));
                        fprintf('%s%s', repmat('', 1, length(progressStr)), newMsg);
                        progressStr = newMsg;
                    end
                    %% Get the image and clean it.
                    try
                        frame = read(v, frameNum);
                    catch
                        continue;
                    end
                    if size(frame, 3) == 3, frame = rgb2gray(frame); end
                    % Detect pupil by thresholding and filling.
                    binaryImage  = (frame > thisParameters.Threshold) & eyeMask;
                    cleanImage   = bwareaopen(imopen(imfill(binaryImage, 'holes'), se), round(thisParameters.Area * obj.MinAreaFrac));


                    %% Find an ellipes or
                    eyePixels  = double(frame(eyeMask));
                    meanEye    = mean(eyePixels);
                    stats_ellipse = regionprops(cleanImage, frame, 'Area', 'Centroid', 'MajorAxisLength', 'MinorAxisLength', 'Eccentricity', 'Orientation', 'MeanIntensity');
                    if ~isempty(stats_ellipse)
                        bestIdx = sbx.PupilTracker.bestRegion(stats_ellipse, thisParameters.StartX, thisParameters.StartY, thisParameters.Area);
                        bestEllipse = stats_ellipse(bestIdx);

                        X(iFrame) = bestEllipse.Centroid(1);
                        Y(iFrame) = bestEllipse.Centroid(2);
                        PupilArea(iFrame) = bestEllipse.Area;

                        MajorAxis(iFrame) = bestEllipse.MajorAxisLength;
                        MinorAxis(iFrame) = bestEllipse.MinorAxisLength;
                        Eccentricity(iFrame) = bestEllipse.Eccentricity;
                        Orientation(iFrame) = bestEllipse.Orientation;

                        ellipseArea = pi * (bestEllipse.MajorAxisLength / 2) * (bestEllipse.MinorAxisLength / 2);
                        if ellipseArea > 0
                            FitQuality(iFrame) = bestEllipse.Area / ellipseArea;
                        end
                        if meanEye > 0
                            EllipseIR(iFrame) = bestEllipse.MeanIntensity / meanEye;
                        end
                    end

                    stats_blob = regionprops(cleanImage, frame, 'Area', 'Centroid', 'BoundingBox', 'MeanIntensity');
                    if ~isempty(stats_blob)
                        bestBlobIdx = sbx.PupilTracker.bestRegion(stats_blob, thisParameters.StartX, thisParameters.StartY, thisParameters.Area);
                        bestBlob = stats_blob(bestBlobIdx);
                        BoundingBox(iFrame,:) = bestBlob.BoundingBox;
                        if meanEye > 0
                            BlobIR(iFrame) = bestBlob.MeanIntensity / meanEye;
                        end
                    end

                    %% Visualise current frame
                    if pv.visualize && ~isempty(obj.Figure) && isgraphics(obj.Figure)
                        ax_v = getappdata(obj.Figure, 'previewAxesHandle');
                        if isempty(ax_v) || ~isgraphics(ax_v)
                            ax_v = axes(obj.Figure, 'Position', [0 0.11 1 0.84], 'Tag', 'previewAxes');
                            setappdata(obj.Figure, 'previewAxesHandle', ax_v);
                        end
                        hold(ax_v, 'off');  % ensure imshow clears the axes on every frame
                        imshow(frame, 'Parent', ax_v, 'initialMagnification', 'fit');
                        hold(ax_v, 'on');
                        plot(ax_v, xy_eye_v(1,:) + thisParameters.EyeX, xy_eye_v(2,:) + thisParameters.EyeY, 'b--', 'LineWidth', 1);
                        if ~isnan(X(iFrame))
                            plot(ax_v, X(iFrame), Y(iFrame), 'g+', 'MarkerSize', 8, 'LineWidth', 2);
                            if ~isnan(MajorAxis(iFrame))
                                phi_p = linspace(0, 2*pi, 50);
                                ap = MajorAxis(iFrame)/2;  bp = MinorAxis(iFrame)/2;
                                tp = pi * Orientation(iFrame) / 180;
                                Rp = [cos(tp) sin(tp); -sin(tp) cos(tp)];
                                xyp = Rp * [ap*cos(phi_p); bp*sin(phi_p)];
                                plot(ax_v, xyp(1,:)+X(iFrame), xyp(2,:)+Y(iFrame), 'r-', 'LineWidth', 2);
                            end
                        end
                        if ~isnan(BoundingBox(iFrame,1))
                            rectangle(ax_v, 'Position', BoundingBox(iFrame,:), 'EdgeColor', 'g', 'LineWidth', 2);
                        end
                        hold(ax_v, 'off');
                        if isnan(X(iFrame))
                            title(ax_v, sprintf('Frame %d | BLINK', frameNum), 'Color', 'r');
                        else
                            title(ax_v, sprintf('Frame %d/%d | Area: %d | Fit: %.2f | EIR: %.2f | BIR: %.2f', ...
                                frameNum, lastFrame, PupilArea(iFrame), FitQuality(iFrame), EllipseIR(iFrame), BlobIR(iFrame)), ...
                                'Interpreter', 'none');
                        end
                        % Update scrubber position and label
                        jSB = getappdata(obj.Figure, 'previewJScrollBar');
                        if ~isempty(jSB)
                            setappdata(obj.Figure, 'previewIgnoreSeekCallback', true);
                            jSB.setValue(max(1, min(frameNum, jSB.getMaximum() - jSB.getVisibleAmount())));
                            setappdata(obj.Figure, 'previewIgnoreSeekCallback', false);
                        end
                        lbl = findobj(obj.Figure, 'Tag', 'previewFrameLabel');
                        if ~isempty(lbl), lbl(1).String = sprintf('%d / %d', frameNum, lastFrame); end
                        setappdata(obj.Figure, 'previewLastRenderedFrame', frameNum);
                        drawnow limitrate;
                    end

                    frameNum = frameNum + pv.frameStep;
                end  % while frameNum <= lastFrame
                nFrames = iFrame;
                fprintf('\n');

                if figClosed, continue; end  % vIdx — skip saving, go to next file

                % Store results in the map
                idx = 1:nFrames;
                resultTable = table(Frames(idx), X(idx), Y(idx), PupilArea(idx), ...
                    MajorAxis(idx), MinorAxis(idx), Eccentricity(idx), Orientation(idx), ...
                    BoundingBox(idx,:),...
                    ThresholdUsed(idx), FitQuality(idx), EllipseIR(idx), BlobIR(idx));
                resultTable.Properties.VariableNames = {'Frame', 'X', 'Y', 'Area', ...
                    'MajorAxis', 'MinorAxis', 'Eccentricity', 'Orientation', ...
                    'BBox','Threshold','FitQuality', 'EllipseIR', 'BlobIR'};
                obj.Results(videoFile) = resultTable;
            end

            if pv.save, obj.save(); end
            if pv.visualize && ~isempty(obj.Figure) && isgraphics(obj.Figure)
                % Tracking finished — start a MATLAB timer to service the scrubber.
                % Java-EDT callbacks (AdjustmentValueChanged) cannot reliably execute
                % MATLAB graphics while the figure is blocked in uiwait.  Timers run on
                % the MATLAB thread and work correctly during uiwait.
                if exist('videoFile','var') && exist('resultTable','var')
                    setappdata(obj.Figure, 'previewVideoFile',   videoFile);
                    setappdata(obj.Figure, 'previewResultTable', resultTable);
                    setappdata(obj.Figure, 'previewParams',      thisParameters);
                    setappdata(obj.Figure, 'previewLastFrame',   lastFrame);
                    setappdata(obj.Figure, 'previewEyeOutline',  xy_eye_v);
                    setappdata(obj.Figure, 'previewMode',        'review');
                    setappdata(obj.Figure, 'previewSeekFrame',   0);
                end
                reviewTimer = timer('ExecutionMode', 'fixedRate', 'Period', 0.15, ...
                    'BusyMode', 'drop', ...
                    'TimerFcn', @(~,~) sbx.PupilTracker.renderStoredFrame(obj.Figure));
                start(reviewTimer);
                uicontrol(obj.Figure, 'Style', 'pushbutton', 'String', 'Close', ...
                    'Units', 'normalized', 'Position', [0.90 0.01 0.09 0.10], ...
                    'FontSize', 10, 'FontWeight', 'bold', ...
                    'Callback', @(~,~) sbx.PupilTracker.stopAndResume(obj.Figure, reviewTimer));
                ax_v = getappdata(obj.Figure, 'previewAxesHandle');
                if ~isempty(ax_v) && isgraphics(ax_v)
                    title(ax_v, 'Tracking complete — scrub to review frames', ...
                        'Color', [0.2 0.6 0.2], 'FontSize', 11);
                end
                drawnow;
                uiwait(obj.Figure);
                if isvalid(reviewTimer), stop(reviewTimer); delete(reviewTimer); end
                if isgraphics(obj.Figure), close(obj.Figure); end
            end
            disp('Batch processing complete!');
        end

        function save(obj)
            % SAVE  Write all tracking results to _pupil.tsv files next to each video.
            %
            %   obj.save()
            %
            % For each entry in obj.Results, writes a tab-separated file whose
            % path is derived by replacing '.mj2' with '_pupil.tsv' in the key.
            %
            % See also  sbx.PupilTracker.track, sbx.PupilTracker.read
            videoFiles = keys(obj.Results);
            for i = 1:numel(videoFiles)
                videoFile = videoFiles{i};
                outputFile = strrep(videoFile, '.mj2', '_pupil.tsv');
                writetable(obj.Results(videoFile), outputFile, 'FileType', 'text', 'Delimiter', '\t');
                fprintf('    -> Saved data to: %s\n', outputFile);
            end
        end

        function read(obj)
            % READ - Load existing _pupil.tsv files into obj.Results.
            videoFiles = keys(obj.Parameters);
            numVideos = numel(videoFiles);
            loaded = 0;
            for vIdx = 1:numVideos
                tsvFile = strrep(videoFiles{vIdx}, '.mj2', '_pupil.tsv');
                if exist(tsvFile, 'file')
                    obj.Results(videoFiles{vIdx}) = readtable(tsvFile, 'FileType', 'text', 'Delimiter', '\t');
                    fprintf('  -> Loaded: %s\n', tsvFile);
                    loaded = loaded + 1;
                else
                    fprintf('  -> Not found, skipping: %s\n', tsvFile);
                end
            end
            fprintf('Done. %d/%d result(s) loaded.\n', loaded, numVideos);
        end

        function plot(obj, columns, pv)
            % PLOT - Plot one or more result columns for each tracked video.
            %   obj.plot("Area")
            %   obj.plot(["Area","FitQuality","EllipseIR","BlobIR"])
            %   obj.plot(["Area","FitQuality"], normalize=false)
            arguments
                obj
                columns (1,:) string {mustBeMember(columns, ["X","Y","Area","MajorAxis","MinorAxis", ...
                    "Eccentricity","Orientation","BBox",...
                    "Threshold","FitQuality","EllipseIR","BlobIR"])} = "Area"
                pv.normalize (1,1) logical = true   % Min-max normalize to [0,1] when multiple columns are plotted
                pv.Parent = []
            end
            doNormalize = pv.normalize && numel(columns) > 1;
            resultFiles = keys(obj.Results);
            if isempty(resultFiles)
                disp('No results to plot. Run track() or read() first.');
                return;
            end
            if isempty(pv.Parent)
                fig = figure('Name', 'PupilTracker Results');
            else
                fig = pv.Parent;
            end
            tiledlayout(fig, numel(resultFiles), 1, 'TileSpacing', 'compact', 'Padding', 'compact');
            for i = 1:numel(resultFiles)
                t = obj.Results(resultFiles{i});
                ax = nexttile;
                hold(ax, 'on');
                for c = 1:numel(columns)
                    col = columns(c);
                    if ismember(col, t.Properties.VariableNames)
                        vals = t.(col);
                        if doNormalize
                            lo = min(vals, [], 'omitnan');
                            hi = max(vals, [], 'omitnan');
                            if hi > lo
                                vals = (vals - lo) ./ (hi - lo);
                            else
                                vals = zeros(size(vals));
                            end
                        end
                        plot(ax, 1:size(vals,1), vals, 'DisplayName', col);
                    else
                        warning('Column "%s" not found in Results for "%s".', col, resultFiles{i});
                    end
                end
                hold(ax, 'off');
                [~, fname] = fileparts(resultFiles{i});
                title(ax, fname, 'Interpreter', 'none');
                xlabel(ax, 'Frame');
                if numel(columns) > 1
                    legend(ax, 'Location', 'best');
                    if doNormalize
                        ylabel(ax, 'Normalized value [0–1]');
                    end
                else
                    ylabel(ax, columns(1));
                end
            end
        end
    end

    methods (Static)
        function pt = fromFolder(rootFolder)
            % FROMFOLDER - Recursively find all *_eye.mj2 files under rootFolder
            % that do not yet have a corresponding *_eye_pupil.json parameter file,
            % and return a PupilTracker configured for those files.
            %
            %   pt = sbx.PupilTracker.fromFolder("D:\data\subject01")
            arguments
                rootFolder (1,1) string
            end
            allMj2  = dir(fullfile(rootFolder, '**', '*_eye.mj2'));
            pending = {};
            for k = 1:numel(allMj2)
                mj2Path  = fullfile(allMj2(k).folder, allMj2(k).name);
                jsonPath = strrep(mj2Path, '.mj2', '_pupil.json');
                if ~isfile(jsonPath)
                    pending{end+1} = mj2Path; %#ok<AGROW>
                end
            end
            if isempty(pending)
                fprintf('No untracked *_eye.mj2 files found under %s\n', rootFolder);
                pt = [];
                return;
            end
            fprintf('Found %d untracked file(s):\n', numel(pending));
            for k = 1:numel(pending)
                fprintf('  %s\n', pending{k});
            end
            % Construct with the first file (constructor takes a scalar string),
            % then add the remaining files directly to the Parameters map.
            pt = sbx.PupilTracker(string(pending{1}));
            if numel(pending) > 1
                emptyRow = table(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
                    'VariableNames', {'Threshold', 'Floor', 'StartX', 'StartY', 'Area', ...
                    'EyeX', 'EyeY', 'EyeMajor', 'EyeMinor', 'EyeOrientation'});
                for k = 2:numel(pending)
                    pt.Parameters(pending{k}) = emptyRow;
                end
            end
        end
    end

    methods (Access = private, Static)

        function resumeWithFlag(fig, flag, val)
            setappdata(fig, flag, val);
            uiresume(fig);
        end

        function restartInitFile(fig)
            % RESTARTINITFILE - Restart button callback: sets the restart flag
            % and resumes uiwait so the current phase unblocks immediately.
            % Disable the button immediately to prevent a rapid double-click
            % from queuing a second uiresume that would bypass Phase 1.
            if isgraphics(fig)
                btn = findobj(fig, 'Style', 'pushbutton', 'String', 'Restart');
                if ~isempty(btn), set(btn, 'Enable', 'off'); end
                setappdata(fig, 'restartFile', true);
                uiresume(fig);
            end
        end

        function roiConfirmIfDouble(fig, evt)
            % Confirm (advance to next step) when the user double-clicks the ROI.
            % images.roi.ROIClickedEventData uses 'double' for double-clicks
            % (distinct from the figure SelectionType which uses 'open').
            if strcmpi(evt.SelectionType, 'double') && isgraphics(fig)
                uiresume(fig);
            end
        end

        function resampleSingleFrame(v, ax)
            % RESAMPLESINGLEFRAME - Button callback for Phase 1:
            % reads one new random frame and swaps the image CData in-place
            % so any drawn ROI is preserved.
            newFrame = sbx.PupilTracker.getFramePool(v,1);
            hImg = findobj(ax, 'Type', 'image');
            if ~isempty(hImg)
                hImg(1).CData = newFrame;
                drawnow;
            end
        end

        function resampleAverageFrame(v, ax, n)
            % RESAMPLEAVERAGEFRAME - Button callback for Phase 1 (average mode):
            % re-samples n random frames, averages them, and swaps CData in-place.
            pool = sbx.PupilTracker.getFramePool(v, n);
            brightnessPerFrame = squeeze(mean(pool,[1 2]));
            keep = brightnessPerFrame > prctile(brightnessPerFrame,50);
            avgFrame = uint8(mean(double(pool(:,:,keep)), 3));
            hImg = findobj(ax, 'Type', 'image');
            if ~isempty(hImg)
                hImg(1).CData = avgFrame;
                drawnow;
            end
        end

        function idx = bestRegion(stats, refX, refY, targetArea)
            % BESTREGION - Select the region with the best composite rank.
            % Regions are ranked on three criteria (lower rank = better):
            %   1) Distance of centroid from (refX, refY)        -> ascending
            %   2) Absolute deviation of area from targetArea     -> ascending
            %   3) Mean intensity                                  -> descending
            % The region with the lowest summed rank is returned.
            n = numel(stats);
            if n == 1, idx = 1; return; end
            centroids = vertcat(stats.Centroid);
            dists  = sqrt((centroids(:,1) - refX).^2 + (centroids(:,2) - refY).^2);
            areas  = [stats.Area]';
            intens = [stats.MeanIntensity]';
            [~, o] = sort(dists,                'ascend');  rankPos  = zeros(n,1); rankPos(o)  = 1:n;
            [~, o] = sort(abs(areas-targetArea),'ascend');  rankArea = zeros(n,1); rankArea(o) = 1:n;
            [~, o] = sort(intens,               'descend'); rankInt  = zeros(n,1); rankInt(o)  = 1:n;
            [~, idx] = min(rankPos + rankArea + rankInt);
        end


        function seekToFrame(fig, jSB)
            % SEEKTOFRAME - Scrubber callback: store the requested frame so the
            % tracking loop picks it up at the next iteration.
            if ~isgraphics(fig), return; end
            if isappdata(fig, 'previewIgnoreSeekCallback') && getappdata(fig, 'previewIgnoreSeekCallback')
                return;
            end
            val = max(jSB.getMinimum(), ...
                min(jSB.getValue(), jSB.getMaximum() - jSB.getVisibleAmount()));
            setappdata(fig, 'previewSeekFrame', val);
        end

        function renderStoredFrame(fig)
            % RENDERSTOREDFRAME - Timer callback used after tracking is complete.
            % Reads the current scrubber position and renders that frame with
            % stored results once the file has finished processing.
            if ~isgraphics(fig), return; end
            if ~isappdata(fig, 'previewMode') || ~strcmp(getappdata(fig, 'previewMode'), 'review')
                return;
            end

            jSB = getappdata(fig, 'previewJScrollBar');
            if isempty(jSB), return; end
            frameNum = max(jSB.getMinimum(), ...
                min(jSB.getValue(), jSB.getMaximum() - jSB.getVisibleAmount()));
            lastRendered = getappdata(fig, 'previewLastRenderedFrame');
            if ~isempty(lastRendered) && isequal(lastRendered, frameNum)
                return;
            end

            videoFile      = getappdata(fig, 'previewVideoFile');
            resultTable    = getappdata(fig, 'previewResultTable');
            thisParameters = getappdata(fig, 'previewParams');
            xy_eye_v       = getappdata(fig, 'previewEyeOutline');
            lastFrame      = getappdata(fig, 'previewLastFrame');
            if isempty(videoFile) || isempty(resultTable), return; end

            try
                vr    = VideoReader(char(videoFile));
                frame = read(vr, frameNum);
            catch
                return;
            end
            if size(frame,3) == 3, frame = rgb2gray(frame); end

            % Find the result row whose Frame value is closest to frameNum
            [~, rowIdx] = min(abs(resultTable.Frame - frameNum));
            X_val     = resultTable.X(rowIdx);
            Y_val     = resultTable.Y(rowIdx);
            areaVal   = resultTable.Area(rowIdx);
            majAx     = resultTable.MajorAxis(rowIdx);
            minAx     = resultTable.MinorAxis(rowIdx);
            orientVal = resultTable.Orientation(rowIdx);
            bboxVal   = resultTable.BBox(rowIdx,:);
            fitQ      = resultTable.FitQuality(rowIdx);
            eirVal    = resultTable.EllipseIR(rowIdx);
            birVal    = resultTable.BlobIR(rowIdx);

            ax_v = getappdata(fig, 'previewAxesHandle');
            if isempty(ax_v) || ~isgraphics(ax_v)
                ax_v = axes(fig, 'Position', [0 0.11 1 0.84], 'Tag', 'previewAxes');
                setappdata(fig, 'previewAxesHandle', ax_v);
            end

            hold(ax_v, 'off');  % ensure imshow clears the axes
            imshow(frame, 'Parent', ax_v, 'InitialMagnification', 'fit');
            hold(ax_v, 'on');
            plot(ax_v, xy_eye_v(1,:) + thisParameters.EyeX, xy_eye_v(2,:) + thisParameters.EyeY, ...
                'b--', 'LineWidth', 1);
            if ~isnan(X_val)
                plot(ax_v, X_val, Y_val, 'g+', 'MarkerSize', 8, 'LineWidth', 2);
                if ~isnan(majAx)
                    phi_p = linspace(0, 2*pi, 50);
                    ap = majAx/2;  bp = minAx/2;
                    tp = pi * orientVal / 180;
                    Rp = [cos(tp) sin(tp); -sin(tp) cos(tp)];
                    xyp = Rp * [ap*cos(phi_p); bp*sin(phi_p)];
                    plot(ax_v, xyp(1,:)+X_val, xyp(2,:)+Y_val, 'r-', 'LineWidth', 2);
                end
            end
            if ~isnan(bboxVal(1))
                rectangle(ax_v, 'Position', bboxVal, 'EdgeColor', 'g', 'LineWidth', 2);
            end
            if isnan(X_val)
                title(ax_v, sprintf('Frame %d | BLINK', frameNum), 'Color', 'r');
            else
                title(ax_v, sprintf('Frame %d/%d | Area: %d | Fit: %.2f | EIR: %.2f | BIR: %.2f', ...
                    frameNum, lastFrame, areaVal, fitQ, eirVal, birVal), 'Interpreter', 'none');
            end
            lbl = findobj(fig, 'Tag', 'previewFrameLabel');
            if ~isempty(lbl), lbl(1).String = sprintf('%d / %d', frameNum, lastFrame); end
            setappdata(fig, 'previewLastRenderedFrame', frameNum);
            drawnow;
        end

        function stopAndResume(fig, t)
            % STOPANDRESUME - Close callback for the post-tracking review window.
            % Stops the review timer and resumes uiwait so track() can return.
            if isvalid(t), stop(t); delete(t); end
            if isgraphics(fig), uiresume(fig); end
        end

        function img = getFramePool(v, n)
            % GETFRAMEPOOL - Sample n random frames from a VideoReader and
            % return them as an H x W x n array.
            totalFrames = round(v.Duration * v.FrameRate);
            n = min(n, totalFrames);
            frameIdx = sort(randperm(totalFrames, n));
            fprintf('Sampling %d frames...\n', n);
            img = zeros(v.Height, v.Width, n, sprintf('uint%d', v.BitsPerPixel));
            for ix = 1:n
                try
                    img(:,:,ix) = read(v, frameIdx(ix));
                catch
                    % Near end of video — leave as zeros
                end
            end
        end

    end
end
