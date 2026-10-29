%{
# Properties of a ROI in a session, based on a Preprocessed set.
-> sbx.Preprocessed
roi : smallint    #  id within this segmentation
---
plane = 0: smallint #  plane
pcell = 0 : float # Probability that the ROI is a neuron
x        : Decimal(4,0) # x-Pixel location
y        : Decimal(4,0) # y-Pixel location
radius   : Decimal(4,1) # Radius of cell in micron
aspect   : Decimal(4,1)  # Aspect raio of major/minor axis of a 2D Gaussian fit.
compact  : Decimal(4,2)  # How compact the ROI is ( 1 is a disk, >1 means less compact)
%}
%
% Note that ROI numbering continues across planes. E.g., with 10 ROIs in
% plane0 and 20 in plane1, roi will run from 1 to 30, with 20:30 in plane1.
% The ROI numbers have to be unique to allow joins with the CChannel table
% that has a single channel number that identifies a time series.
classdef PreprocessedRoi < dj.Part
    properties (SetAccess = protected)
        master = sbx.Preprocessed
    end

    methods  (Access=public)
        function nwbRoot = nwb(tbl,nwbRoot,pv)
            % Add to NWB root
            % Read the ROI masks from stat file
            % Read the associated fluorescence data from CChannel tbl and
            % store

            prep = sbx.Preprocessed & tbl; %
            assert(count(prep)==1,"More than one preprocessed set.NIY");
            stat =  prep.stat; % Reads the stat.npy file
            assert(~isempty(stat),'Stat file not found.');
            nrRoi = count(tbl);

            % Validate that all ROIs have pixel data
            npixArray = [stat.npix];
            if any(npixArray == 0)
                warning('Some ROIs have no pixels, skipping NWB export');
                return;
            end

            nrPix = sum(npixArray);
            % Create proper pixel mask format for NWB compound type
            % NWB expects a table-like structure with columns [x, y, weight]
            allY = [stat.ypix];
            allX = [stat.xpix];

            % Ensure we have the expected number of pixels
            if length(allY) ~= nrPix || length(allX) ~= nrPix
                error('Mismatch in pixel count: expected %d, got Y:%d, X:%d', nrPix, length(allY), length(allX));
            end

            % Create the pixel mask as a proper table format
            % Note: NWB pixel masks are typically stored as [x, y, weight] not [y, x, weight]
            pixelMaskData = table(...
                uint32(allX(:)), ...
                uint32(allY(:)), ...
                single(ones(nrPix,1)), ...
                'VariableNames', {'x', 'y', 'weight'});

            % Index where each ROI ends
            maskIndex = uint64(cumsum(npixArray));

            % Create VectorData and VectorIndex for the pixel masks
            pixelMask = types.hdmf_common.VectorData(...
                'data', pixelMaskData, ...
                'description', 'Pixel masks for ROIs');
            pixelMaskIndex = types.hdmf_common.VectorIndex(...
                'target', types.untyped.ObjectView(pixelMask), ...
                'data', maskIndex);


            % Add to a plane segmentation
            imgPlane = get(nwbRoot.general_optophysiology,'imaging_plane');
            planeSegmentation = types.core.PlaneSegmentation(...
                'colnames', {'pixel_mask'}, ...
                'description', 'Segmented by Suite2p', ...
                'imaging_plane', types.untyped.SoftLink(imgPlane), ...
                'pixel_mask', pixelMask, ...
                'pixel_mask_index', pixelMaskIndex);
            imgSegmentation  = types.core.ImageSegmentation();
            imgSegmentation.planesegmentation.set('PlaneSegmentation',planeSegmentation);

            ophys = types.core.ProcessingModule('description','Optophysiology');
            ophys.nwbdatainterface.set('ImageSegmentation',imgSegmentation);
            nwbRoot.processing.set('ophys',ophys);

            % Collect the signals from the nsCChannel table
            % Here we have to restrict by the pv.experiment to pull only
            % the data for this experiment and not for all experiments in
            % the session.
            ctag = 'fluorescence';
            cChannelTbl = (ns.CChannel & struct('ctag',ctag) & (ns.Experiment & pv.experiment)) & proj(tbl,'subject','session_date','roi->channel');
            nrRoiToExport = count(cChannelTbl);
            if  nrRoiToExport ~=nrRoi
                fprintf("Exporting %s for %d out of %d ROIs \n",ctag, nrRoiToExport ,nrRoi);                
            end
            [signal,roiNr] = fetchn(cChannelTbl,'signal','channel');
            
            signal = cat(2,signal{:})'; %[nrRoi nrFrames]
            signal   = types.untyped.DataPipe('data',signal,'ChunkSize',[nrRoiToExport  1]);
            cTpl = fetch(ns.C&cChannelTbl,'time','ctag');
            assert(numel(cTpl.time)==3,'Time not regularly sampled?');
            nstime =fetch1(ns.PluginParameter  & (ns.Experiment &pv.experiment) & 'plugin_name="cic"' & 'property_name="trial"' ,'property_nstime');
            timeZero  = nstime(1);% Start of first trial =0
            startTime = (cTpl.time(1)-timeZero)/1000; % Align and convert to seconds
            sampleRate = 1000*(cTpl.time(3)-1)./(cTpl.time(2)-cTpl.time(1)); % Hertz
            roiTableRegion = types.hdmf_common.DynamicTableRegion( ...
                'table',types.untyped.ObjectView(planeSegmentation), ...
                'description','all rois', ...
                'data',int64(roiNr-1));
            roiResponseSeries = types.core.RoiResponseSeries  ( ...
                'rois',roiTableRegion, ...
                'data',signal, ...
                'data_unit','spk/s', ...
                'starting_time',startTime, ...
                'starting_time_rate',sampleRate, ...
                'description',cTpl.ctag);

            fluorescence = types.core.Fluorescence();
            fluorescence.roiresponseseries.set('RoiResponseSeries',roiResponseSeries);
            ophys.nwbdatainterface.set('Fluorescence',fluorescence);
        end
        function [hData] = plotSpatial(roi,pv)
            % Show properties of ROIs in a spatial layout matching that
            % used in suite2p.
            %
            % INPUT
            % roi  = (subset of) sbx.Roi
            %
            % Parameter/Value pairs:
            % sz = Size to use for each ROI. By default, the estimated physical size
            %       of the roi is used. But by passing some other property (e.g. tuning per ROI),
            %       this can be visualized. Maximum size is 500 points (surface area).
            % szLabel  - Label to use in the datatip.
            % color = Color to use for each ROI. By default this represents the compactness.
            % colorLabel = Label to use in the datatip.
            % clim - Color limits to use
            % showImg  - Set to true to show the mean Ca image as a
            %               background [true]
            % colormap - Colormap to use  [hot]
            % pix - Set to true to use the actual ROI pixels (instead of a
            % circle)  [false]
            % alpha - Transparency [0.5]
            % alphaThreshold - 0
            % randColor - Assign each ROI a random color from the colormap [false]
            %
            % OUTPUT
            %    hScatter - Handle to the scatter object
            % EXAMPLE
            %
            % Anatimical properties :
            % Size of the circle represents the size of the ROI, color
            % represents the compactness.
            % plot(sbx.PreprocessedRoi & 'session_date="2024-02-09"');
            % Functional properties for a subset of ROI
            % roi = sbx.PreprocessedRoi & 'pcell>0.75';
            %  Query the CChannel table, for instance spike rate in a
            %  specific experiment:
            % cChannel = ns.CChannel & 'ctag="spikes"' & 'starttime="14:08:06"'
            % Then do an inner join to match roi with channels:
            % spikes = fetchn(cChannel & proj(roi,'subject','session_date','roi->channel'),'signal')
            % meanSpikes = cellfun(@mean,spikes)
            % then plot:
            % plotSpatial(roi,color=meanSpikes,pix=true)
            arguments
                roi (1,1) sbx.PreprocessedRoi
                pv.sz   (1,:)  double  = []
                pv.szLabel (1,1) string = ''
                pv.color   (1,:) double = []
                pv.colorLabel (1,1) string = ''
                pv.clim    (1,:) double  =[];
                pv.showImg (1,1) logical = true;
                pv.colormap (:,3) double = hot;
                pv.pix (1,1) logical  = false
                pv.alpha (1,1) double = 0.5;
                pv.alphaThreshold (1,1) double =0;
                pv.randColor (1,1) logical = false;
            end
            nrRoi = count(roi);
            if nrRoi ==0
                fprintf('No rois, nothing to plot \n');
                return;
            end
            % To visualize an arbitrary and large subset of rois, this
            % function can be called with a large struct array (one per roi),
            % which is potentially inefficient. Silencing DJ warning about this.
            warnState = warning('query');
            warning('off', 'DataJoint:longCondition');

            MAXPOINTS = 300; %
            [x,y,radius,roiNr] = fetchn(roi,'x','y','radius','roi');

            [xScale,yScale] =fetchn(sbx.Preprocessed & roi,'xscale','yscale');
            micPerPix = [xScale yScale];
            micPerPixR = sqrt(sum(micPerPix.^2));


            %% Setup color
            if pv.randColor
                % Use discrete colormap entries so each ROI keeps the same
                % color in both scatter and pixel-overlay modes.
                pv.color = randi(size(pv.colormap,1),nrRoi,1);
                pv.colorLabel = 'Random color';
                if isempty(pv.clim)
                    pv.clim = [1 size(pv.colormap,1)];
                end
            elseif isempty(pv.color)
                % Show the compactness as the color of the cells
                pv.color = fetchn(roi,'compact');
                pv.colorLabel = 'Compactness';
            end
            % Clamp to the limits
            if isempty(pv.clim)
                pv.clim = [min(pv.color) max(pv.color)];
            else
                pv.color = max(min(pv.color,pv.clim(2)),pv.clim(1)); % Clamp between clims
            end
            colormap(pv.colormap)


            %% Show mean image as background
            if pv.showImg || pv.pix
                meanImg = fetch1(sbx.Preprocessed & roi,'img','LIMIT 1');
                RI = imref2d(size(meanImg),micPerPix(1),micPerPix(2));
                cl = [min(meanImg,[],"all") max(meanImg,[],"all")];
                meanImg  = 1+round((size(gray,1)-1)*(meanImg-cl(1))./(cl(2)-cl(1)));
                % Convert to true color to reuse axes for color map
                meanImg = ind2rgb(meanImg,gray);
                if pv.showImg 
                imshow(meanImg,RI);
                hold on
                end
            end

            if pv.pix
                %% Show the data in the ROI pixels
                overlay = zeros(RI.ImageSize);
                alpha   = zeros(RI.ImageSize);
                prep = sbx.Preprocessed & roi;
                stat =  prep.stat; % Reads the stat.npy file
                assert(~isempty(stat),'Stat file not loaded. (is NS_ROOT set correctly?) Cannot use pix=true');

                rCntr  =0;
                % Scale to the full colormap
                if pv.randColor
                    clrIndex = pv.color;
                else
                    scaled  = (pv.color-min(pv.color))./(max(pv.color)-min(pv.color));
                    clrIndex  = 1+round((size(pv.colormap,1)-1)*scaled);
                end
                for r=roiNr'
                    rCntr= rCntr+1;
                    x =stat(r).xpix;
                    y = stat(r).ypix;
                    ix =sub2ind(size(overlay),y,x);
                    overlay(ix)= clrIndex(rCntr);
                    alpha(ix) = pv.alpha*(abs(pv.color(rCntr))>pv.alphaThreshold);
                end
                overlay = ind2rgb(overlay,pv.colormap); % Convert to true color using the colormap
                hData = imshow(overlay,RI);
                set(hData,'AlphaData',alpha)
            else
                %% Show the data as circles
                % Setup size
                if isempty(pv.sz)
                    % Use the physical size as the size of the cells
                    pv.sz =min(pi*(radius*micPerPixR).^2,MAXPOINTS);
                    pv.szLabel = 'size (\mum^2)';
                end
                zeroSize = pv.sz==0;
                if any(zeroSize)
                    pv.sz(zeroSize) = eps;
                end

                % Show the data symbolically on top of the image. Match
                % orientation to suite2p
                hData = scatter((x-1)*micPerPix(1),(y-1)*micPerPix(2),pv.sz(:), pv.color(:),'filled','MarkerFaceAlpha',pv.alpha);%.*(abs(pv.color(:))>pv.alphaThreshold));
                set(gca,'Color','none','XDir','normal','YDir','reverse','colormap',pv.colormap);

                % Data tips
                if ~isempty(pv.szLabel)
                    hData.DataTipTemplate.DataTipRows(3).Label = pv.szLabel;
                end
                if ~isempty(pv.colorLabel)
                    hData.DataTipTemplate.DataTipRows(4).Label = pv.colorLabel;
                end
            end
            if ~isempty(pv.clim)
                clim(pv.clim)
            end

            %% Add Markup
            xlabel 'X (\mum)'
            ylabel 'Y (\mum)'
            h = colorbar;
            ylabel(h,pv.colorLabel);


            % Restore warning state
            warning(warnState);

        end

        function [movieData] = plotMovie(roi,timeCourse,pv)
            % Show a movie of ROI values in the spatial layout.
            %
            % INPUT
            % roi        = (subset of) sbx.PreprocessedRoi
            % timeCourse = vector of values over time, applied to every ROI,
            %              or an [nRoi x nTime] matrix with one row per ROI.
            %              An [nTime x nRoi] matrix is also accepted.
            %
            % Parameter/Value pairs:
            % clim            - Color limits shared by all frames [min max]
            % colormap        - Colormap to use [hot]
            % showImg         - Show the mean Ca image as background [true]
            % alpha           - ROI transparency [0.5]
            % alphaThreshold  - Values with smaller absolute magnitude are
            %                   made transparent [0]
            % frameRate       - Playback rate in frames per second [10]
            % loopCount       - Number of playback repetitions [1]
            %
            % OUTPUT
            % movieData - MATLAB movie frame structure returned by getframe.
            arguments
                roi (1,1) sbx.PreprocessedRoi
                timeCourse double
                pv.clim (1,:) double = []
                pv.colormap (:,3) double = hot
                pv.showImg (1,1) logical = true
                pv.alpha (1,1) double = 0.5
                pv.alphaThreshold (1,1) double = 0
                pv.frameRate (1,1) double = 30
                pv.loopCount (1,1) double = 1
            end

            nrRoi = count(roi);
            if nrRoi == 0
                error('Cannot create a movie for an empty ROI relation.');
            end
            if isempty(timeCourse)
                error('timeCourse must contain at least one value.');
            end

            % Normalize the time-course to [nRoi x nTime]. A vector is
            % shared by all ROIs; matrices may be supplied in either common
            % orientation.
            if isvector(timeCourse)
                values = repmat(timeCourse(:)',nrRoi,1);
            elseif size(timeCourse,1) == nrRoi
                values = timeCourse;
            elseif size(timeCourse,2) == nrRoi
                values = timeCourse';
            else
                error('timeCourse must be a vector or have one dimension equal to the number of ROIs (%d).',nrRoi);
            end
            nrFrames = size(values,2);

            [x,y,radius,roiNr] = fetchn(roi,'x','y','radius','roi'); %#ok<ASGLU>
            [xScale,yScale] = fetchn(sbx.Preprocessed & roi,'xscale','yscale');
            micPerPix = [xScale yScale];

            prep = sbx.Preprocessed & roi;
            stat = prep.stat;
            assert(~isempty(stat),'Stat file not loaded. (is NS_ROOT set correctly?)');
            meanImg = fetch1(prep,'img','LIMIT 1');
            RI = imref2d(size(meanImg),micPerPix(1),micPerPix(2));

            if isempty(pv.clim)
                pv.clim = [min(values,[],'all') max(values,[],'all')];
            end
            if pv.clim(1) == pv.clim(2)
                pv.clim = pv.clim + [-0.5 0.5];
            end

            if pv.showImg
                imageClim = [min(meanImg,[],'all') max(meanImg,[],'all')];
                imageIndex = 1+round((size(gray,1)-1) * ...
                    (meanImg-imageClim(1)) ./ (imageClim(2)-imageClim(1)));
                background = ind2rgb(imageIndex,gray);
                imshow(background,RI);
                hold on
            else
                imshow(zeros([size(meanImg) 3]),RI);
                hold on
            end

            overlay = zeros([size(meanImg) 3]);
            alphaData = zeros(size(meanImg));
            hData = imshow(overlay,RI);
            set(hData,'AlphaData',alphaData);

            movieData(1,nrFrames) = getframe(gca);
            for frameNr = 1:nrFrames
                colorIndex = 1+round((size(pv.colormap,1)-1) * ...
                    (values(:,frameNr)-pv.clim(1)) ./ diff(pv.clim));
                colorIndex = max(1,min(size(pv.colormap,1),colorIndex));

                overlay = zeros([size(meanImg) 3]);
                alphaData = zeros(size(meanImg));
                for roiCntr = 1:nrRoi
                    thisStat = stat(roiNr(roiCntr));
                    pixelIndex = sub2ind(size(meanImg),thisStat.ypix,thisStat.xpix);
                    for colorChannel = 1:3
                        channel = overlay(:,:,colorChannel);
                        channel(pixelIndex) = pv.colormap(colorIndex(roiCntr),colorChannel);
                        overlay(:,:,colorChannel) = channel;
                    end
                    alphaData(pixelIndex) = pv.alpha * ...
                        (abs(values(roiCntr,frameNr)) > pv.alphaThreshold);
                end

                set(hData,'CData',overlay,'AlphaData',alphaData);
                title(sprintf('Frame %d/%d',frameNr,nrFrames));
                drawnow;
                movieData(frameNr) = getframe(gca);
            end

            movie(gcf,movieData,pv.loopCount,pv.frameRate);
        end
    end
    methods (Access=?sbx.Preprocessed)
        function makeTuples(tbl,key)
            assert(setupPython,"Could not find a Python installation");  % Make sure we have a python environment

            prep = fetch(sbx.Preprocessed& key,'*');
            micPerPix = sqrt(sum([prep.xscale prep.yscale].^2));
            fldr = unique(fullfile(folder(ns.Experiment & key),prep.folder));
            if ~exist(fldr,"dir")
                error('Preprocessed data folder %s not found',fldr)
            end
           % depth0..n has separate recordings at different depths
           % and plane0..n  has simultaneous recordings at different
           % dephts. If there is only a single depth, the depth0 folder is
           % omitted.
            planes = dir(fullfile(fldr,'plane*')); 
            roisSoFar = 0;

            for pl = 1:numel(planes)
                %% Read npy
                thisFile = fullfile(planes(pl).folder,planes(pl).name,'iscell.npy');
                if ~exist(thisFile,"file")
                    error('File %s does not exist',thisFile);
                end
                iscell = ndarrayToArray(py.numpy.load(thisFile,allow_pickle=true),single=true);

                % We saved the stat.npy as stat.mat in sbx.Preprocessed
                thisFile = fullfile(planes(pl).folder,planes(pl).name,'stat.mat');
                if ~exist(thisFile,"file")
                    error('File %s does not exist',thisFile);
                end
                load(thisFile,'stat');
                stat= [stat{:}];
                med= cat(1,stat.med); %[y x] center pixels per ROI.
                compact = cat(1,stat.compact);
                aspect = cat(1,stat.aspect_ratio);

                radius = cat(1,stat.radius); % Pixels
                radius = radius.*micPerPix;

                %% Make tuples and insert=
                nrROIs= numel(aspect);
                roi = roisSoFar+(1:nrROIs);
                tpl = mergestruct(key, ...
                    struct('roi',num2cell(roi)', ...
                    'plane',pl-1, ...
                    'pcell',num2cell(iscell(:,2)), ...
                    'x',num2cell(med(:,2)), ...    % Dim 2 is the horizontal axis of the image
                    'y',num2cell(med(:,1)), ...    % Dim 1 is the vertical axis of the image
                    'radius',num2cell(radius',1)', ...
                    'compact',num2cell(compact',1)', ...
                    'aspect',num2cell(aspect',1)'));
                insert(tbl,tpl);
                roisSoFar = roisSoFar + nrROIs; % Next plane starts at nrROIs+1
            end
        end
    end


end
