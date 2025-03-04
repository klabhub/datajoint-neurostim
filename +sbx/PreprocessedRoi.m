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
        function o = PreprocessedRoi(varargin)
            assert(setupPython,"Could not find a Python installation");  % Make sure we have a python environment
             o@dj.Part(varargin{:});
        end
        function nwbRoot = nwb(tbl,nwbRoot,pv)
            % Add to NWB root
            % Read the masks
              error('Still to be changed after change to dj.Part table')
%%   NOTES:
% This table is no longer linked to an experiment, but to the session.
% Have to adjust nwbExport in ns.Experiment to include this information.
% and separately export some or all of the C table entries (which were
% previously linked via sbx.Roi)

            prep = sbx.Preprocessed & tbl; %
            assert(count(prep)==1,"More than one preprocessed set.NIY");
            stat =  prep.stat; % Reads the stat.npy file
            assert(~isempty(stat),'Stat file not found.');
            sz = size(fetch1(prep,'img'));           
            nrRoi = count(tbl);
            masks =zeros([fliplr(sz) nrRoi]);
            % Loop over the roi, setting pixels in the mask to 1.
            for r=1:nrRoi
                x =stat(r).xpix;
                y = stat(r).ypix;
                roi = r*ones(1,stat(r).npix);
                ix= sub2ind(size(masks),x,y,roi);
                masks(ix) =1;
            end
            % Add to a plane segmentation
            imgPlane = get(nwbRoot.general_optophysiology,'imaging_plane');
            planeSegmentation = types.core.PlaneSegmentation(...
                'colnames',{'image_mask'}, ...
                'description','Segmented by ', ...
                'imaging_plane',types.untyped.SoftLink(imgPlane),...
                'image_mask',types.hdmf_common.VectorData('data',masks,'description','roi masks'));

            imgSegmentation  = types.core.ImageSegmentation();
            imgSegmentation.planesegmentation.set('PlaneSegmentation',planeSegmentation);

            ophys = types.core.ProcessingModule('description','Optophysiology');
            ophys.nwbdatainterface.set('ImageSegmentation',imgSegmentation);
            nwbRoot.processing.set('ophys',ophys);

            % Collect the signals from the nsCChannel table
            cChannelTbl = (ns.CChannel & 'ctag="fluorescence"' & expt) & roi.proj('subject','session_date','roi->channel');
            
            signal = fetchn(ns.CChannel&tbl,'signal');
            signal = cat(2,signal{:})'; %[nrRoi nrFrames]
            signal   = types.untyped.DataPipe('data',signal,'ChunkSize',[nrRoi 1]);
            cTpl = fetch(ns.C&tbl,'time','ctag');
            assert(numel(cTpl.time)==3,'Time not regularly sampled?');
            nstime =fetch1(ns.PluginParameter  & (ns.Experiment &tbl) & 'plugin_name="cic"' & 'property_name="trial"' ,'property_nstime');
            timeZero  = nstime(1);% Start of first trial =0
            startTime = (cTpl.time(1)-timeZero)/1000; % Align and convert to seconds
            sampleRate = 1000*cTpl.time(3)./(cTpl.time(2)-cTpl.time(1)); % Hertz
            roiTableRegion = types.hdmf_common.DynamicTableRegion( ...
                'table',types.untyped.ObjectView(planeSegmentation), ...
                'description','all rois', ...
                'data',(0:nrRoi-1)');
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
            if isempty(pv.color) 
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
            if pv.showImg
                meanImg = fetch1(sbx.Preprocessed & roi,'img','LIMIT 1');
                RI = imref2d(size(meanImg),micPerPix(1),micPerPix(2));
                cl = [min(meanImg,[],"all") max(meanImg,[],"all")];
                meanImg  = 1+round((size(gray,1)-1)*(meanImg-cl(1))./(cl(2)-cl(1)));
                % Convert to true color to reuse axes for color map
                meanImg = ind2rgb(meanImg,gray);
                imshow(meanImg,RI);
                hold on
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
                scaled  = (pv.color-min(pv.color))./(max(pv.color)-min(pv.color));
                clrIndex  = 1+round((size(pv.colormap,1)-1)*scaled);
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
    end
    methods (Access=?sbx.Preprocessed)
        function makeTuples(tbl,key)
            prep = fetch(sbx.Preprocessed& key,'*');
            micPerPix = sqrt(sum([prep.xscale prep.yscale].^2));           
            fldr = unique(fullfile(folder(ns.Experiment & key),prep.folder));
            if ~exist(fldr,"dir")
                error('Preprocessed data folder %s not found',fldr)
            end
            planes = dir(fullfile(fldr,'plane*'));
            roisSoFar = 0;

            for pl = 1:numel(planes)
                %% Read npy
                thisFile = fullfile(fldr,planes(pl).name,'iscell.npy');
                if ~exist(thisFile,"file")
                    error('File %s does not exist',thisFile);
                end
                iscell = ndarrayToArray(py.numpy.load(thisFile,allow_pickle=true),single=true);

                % We saved the stat.npy as stat.mat in sbx.Preprocessed
                thisFile = fullfile(fldr,planes(pl).name,'stat.mat');
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